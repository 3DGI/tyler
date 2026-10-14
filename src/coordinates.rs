use crate::spatial_structs::Bbox;
use cityjson_lib::ops::Transformer;

#[derive(Clone, Copy, Debug)]
pub struct RootEnuFrame {
    pub source_origin: [f64; 3],
    pub geodetic_origin: [f64; 3],
    pub ecef_origin: [f64; 3],
    pub east: [f64; 3],
    pub north: [f64; 3],
    pub up: [f64; 3],
}

impl RootEnuFrame {
    pub fn from_bbox(source_crs: &str, bbox: &Bbox) -> Result<Self, Box<dyn std::error::Error>> {
        let source_origin = [
            f64::midpoint(bbox[0], bbox[3]),
            f64::midpoint(bbox[1], bbox[4]),
            f64::midpoint(bbox[2], bbox[5]),
        ];
        Self::from_source_origin(source_crs, source_origin)
    }

    pub fn from_source_origin(
        source_crs: &str,
        source_origin: [f64; 3],
    ) -> Result<Self, Box<dyn std::error::Error>> {
        let to_geodetic = cityjson_lib::ops::transformer(source_crs, "EPSG:4979")?;
        let to_ecef = cityjson_lib::ops::transformer(source_crs, "EPSG:4978")?;

        let geodetic_origin = to_geodetic.transform(source_origin)?;
        let ecef_origin = to_ecef.transform(source_origin)?;
        let (east, north, up) = enu_basis(geodetic_origin[0], geodetic_origin[1]);

        Ok(Self {
            source_origin,
            geodetic_origin,
            ecef_origin,
            east,
            north,
            up,
        })
    }

    pub fn transform(&self) -> [f64; 16] {
        [
            self.east[0],
            self.east[1],
            self.east[2],
            0.0,
            self.north[0],
            self.north[1],
            self.north[2],
            0.0,
            self.up[0],
            self.up[1],
            self.up[2],
            0.0,
            self.ecef_origin[0],
            self.ecef_origin[1],
            self.ecef_origin[2],
            1.0,
        ]
    }
}

pub fn horizontal_footprint_width_meters(
    source_crs: &str,
    bbox: &Bbox,
    source_to_ecef: &Transformer,
) -> Result<f64, Box<dyn std::error::Error>> {
    let tile_enu_frame = RootEnuFrame::from_bbox(source_crs, bbox)?;
    let z = f64::midpoint(bbox[2], bbox[5]);
    let corners = [
        [bbox[0], bbox[1], z],
        [bbox[0], bbox[4], z],
        [bbox[3], bbox[1], z],
        [bbox[3], bbox[4], z],
    ];
    let mut min_east = f64::INFINITY;
    let mut min_north = f64::INFINITY;
    let mut max_east = f64::NEG_INFINITY;
    let mut max_north = f64::NEG_INFINITY;

    for source_position in corners {
        let ecef_position = source_to_ecef.transform(source_position)?;
        let delta = [
            ecef_position[0] - tile_enu_frame.ecef_origin[0],
            ecef_position[1] - tile_enu_frame.ecef_origin[1],
            ecef_position[2] - tile_enu_frame.ecef_origin[2],
        ];
        let east = dot_product(delta, tile_enu_frame.east);
        let north = dot_product(delta, tile_enu_frame.north);
        min_east = min_east.min(east);
        min_north = min_north.min(north);
        max_east = max_east.max(east);
        max_north = max_north.max(north);
    }

    let width = (max_east - min_east).max(max_north - min_north);
    if !width.is_finite() || width < 0.0 {
        return Err("could not calculate horizontal footprint width in metres".into());
    }
    Ok(width)
}

fn dot_product(lhs: [f64; 3], rhs: [f64; 3]) -> f64 {
    lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2]
}

fn enu_basis(longitude_degrees: f64, latitude_degrees: f64) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let lon = longitude_degrees.to_radians();
    let lat = latitude_degrees.to_radians();
    let (sin_lon, cos_lon) = lon.sin_cos();
    let (sin_lat, cos_lat) = lat.sin_cos();

    let east = [-sin_lon, cos_lon, 0.0];
    let north = [-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat];
    let up = [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat];

    (east, north, up)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dot(lhs: [f64; 3], rhs: [f64; 3]) -> f64 {
        lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2]
    }

    fn cross(lhs: [f64; 3], rhs: [f64; 3]) -> [f64; 3] {
        [
            lhs[1] * rhs[2] - lhs[2] * rhs[1],
            lhs[2] * rhs[0] - lhs[0] * rhs[2],
            lhs[0] * rhs[1] - lhs[1] * rhs[0],
        ]
    }

    fn length(vector: [f64; 3]) -> f64 {
        dot(vector, vector).sqrt()
    }

    fn assert_approx_eq(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1.0e-12,
            "expected {expected}, got {actual}"
        );
    }

    #[test]
    fn root_enu_basis_is_orthonormal_and_right_handed() {
        let frame = RootEnuFrame::from_source_origin("EPSG:4979", [4.9, 52.37, 12.0])
            .expect("EPSG:4979 origin should produce an ENU frame");

        assert_approx_eq(length(frame.east), 1.0);
        assert_approx_eq(length(frame.north), 1.0);
        assert_approx_eq(length(frame.up), 1.0);

        assert_approx_eq(dot(frame.east, frame.north), 0.0);
        assert_approx_eq(dot(frame.east, frame.up), 0.0);
        assert_approx_eq(dot(frame.north, frame.up), 0.0);

        let east_cross_north = cross(frame.east, frame.north);
        for (actual, expected) in east_cross_north.into_iter().zip(frame.up) {
            assert_approx_eq(actual, expected);
        }
    }

    #[test]
    fn root_enu_transform_translation_is_ecef_origin() {
        let frame = RootEnuFrame::from_source_origin("EPSG:4979", [4.9, 52.37, 12.0])
            .expect("EPSG:4979 origin should produce an ENU frame");
        let transform = frame.transform();

        assert_approx_eq(transform[12], frame.ecef_origin[0]);
        assert_approx_eq(transform[13], frame.ecef_origin[1]);
        assert_approx_eq(transform[14], frame.ecef_origin[2]);
        assert_approx_eq(transform[15], 1.0);
    }
}
