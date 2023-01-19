//! Proj.convert() function adapted from the
//! [proj](https://github.com/georust/proj/blob/main/src/proj.rs) crate to transform xyz
//! coordinates instead of only xy.
//!
//! TODO: include proj license

use libc::c_int;
use libc::{c_char, c_double};
use num_traits::Float;
use proj_sys::{
    proj_area_create, proj_area_set_bbox, proj_context_create, proj_context_errno,
    proj_create_crs_to_crs, proj_destroy, proj_errno_string, proj_normalize_for_visualization,
    proj_trans, PJconsts, PJ_AREA, PJ_CONTEXT, PJ_COORD, PJ_DIRECTION_PJ_FWD, PJ_XYZT,
};
use std::{
    fmt::{self, Debug},
    str,
};

use proj_sys::{proj_errno, proj_errno_reset};

use std::ffi::{CStr, CString, NulError};
use thiserror::Error;

pub trait CoordinateType: Float + Copy + PartialOrd + Debug {}

impl<T: Float + Copy + PartialOrd + Debug> CoordinateType for T {}

/// Called by new_known_crs and proj_known_crs
fn transform_epsg(
    ctx: *mut PJ_CONTEXT,
    from: &str,
    to: &str,
    area: Option<Area>,
) -> Result<Proj, ProjCreateError> {
    let from_c = CString::new(from).map_err(ProjCreateError::ArgumentNulError)?;
    let to_c = CString::new(to).map_err(ProjCreateError::ArgumentNulError)?;
    let proj_area = unsafe { proj_area_create() };
    area_set_bbox(proj_area, area);
    let ptr = result_from_create(ctx, unsafe {
        proj_create_crs_to_crs(ctx, from_c.as_ptr(), to_c.as_ptr(), proj_area)
    })
    .map_err(|e| ProjCreateError::ProjError(e.message(ctx)))?;
    // Normalise input and output order to Lon, Lat / Easting Northing by inserting
    // An axis swap operation if necessary
    let normalised = unsafe {
        let normalised = proj_normalize_for_visualization(ctx, ptr);
        // deallocate stale PJ pointer
        proj_destroy(ptr);
        normalised
    };
    Ok(Proj {
        c_proj: normalised,
        ctx,
        area: Some(proj_area),
    })
}

/// Construct a `Result` from the result of a `proj_create*` call.
fn result_from_create<T>(context: *mut PJ_CONTEXT, ptr: *mut T) -> Result<*mut T, Errno> {
    if ptr.is_null() {
        Err(Errno(unsafe { proj_context_errno(context) }))
    } else {
        Ok(ptr)
    }
}

#[derive(Copy, Clone, Debug)]
pub struct Area {
    pub north: f64,
    pub south: f64,
    pub east: f64,
    pub west: f64,
}

fn area_set_bbox(parea: *mut proj_sys::PJ_AREA, new_area: Option<Area>) {
    // if a bounding box has been passed, modify the proj area object
    if let Some(narea) = new_area {
        unsafe {
            proj_area_set_bbox(parea, narea.west, narea.south, narea.east, narea.north);
        }
    }
}

pub trait Coord<T>
where
    T: CoordinateType,
{
    fn x(&self) -> T;
    fn y(&self) -> T;
    fn z(&self) -> T;
    fn from_xyz(x: T, y: T, z: T) -> Self;
}

impl<T: CoordinateType> Coord<T> for (T, T, T) {
    fn x(&self) -> T {
        self.0
    }
    fn y(&self) -> T {
        self.1
    }
    fn z(&self) -> T {
        self.2
    }
    fn from_xyz(x: T, y: T, z: T) -> Self {
        (x, y, z)
    }
}

pub struct Proj {
    c_proj: *mut PJconsts,
    ctx: *mut PJ_CONTEXT,
    area: Option<*mut PJ_AREA>,
}

impl Proj {
    pub fn new_known_crs(
        from: &str,
        to: &str,
        area: Option<Area>,
    ) -> Result<Proj, ProjCreateError> {
        let ctx = unsafe { proj_context_create() };
        transform_epsg(ctx, from, to, area)
    }

    pub fn convert<C, F>(&self, point: C) -> Result<C, ProjError>
    where
        C: Coord<F>,
        F: CoordinateType,
    {
        let c_x: c_double = point.x().to_f64().ok_or(ProjError::FloatConversion)?;
        let c_y: c_double = point.y().to_f64().ok_or(ProjError::FloatConversion)?;
        let c_z: c_double = point.z().to_f64().ok_or(ProjError::FloatConversion)?;
        let new_x;
        let new_y;
        let new_z;
        let err;

        // This doesn't seem strictly correct, but if we set PJ_XY or PJ_LP here, the
        // other two values remain uninitialized and we can't be sure that libproj
        // doesn't try to read them. proj_trans_generic does the same thing.
        let xyzt = PJ_XYZT {
            x: c_x,
            y: c_y,
            z: c_z,
            t: f64::INFINITY,
        };
        unsafe {
            proj_errno_reset(self.c_proj);
            let trans = proj_trans(self.c_proj, PJ_DIRECTION_PJ_FWD, PJ_COORD { xyzt });
            new_x = trans.xyz.x;
            new_y = trans.xyz.y;
            new_z = trans.xyz.z;
            err = proj_errno(self.c_proj);
        }
        if err == 0 {
            Ok(C::from_xyz(
                F::from(new_x).ok_or(ProjError::FloatConversion)?,
                F::from(new_y).ok_or(ProjError::FloatConversion)?,
                F::from(new_z).ok_or(ProjError::FloatConversion)?,
            ))
        } else {
            Err(ProjError::Conversion(error_message(err)?))
        }
    }
}

/// Errors originating in PROJ which can occur during projection and conversion
#[derive(Error, Debug)]
pub enum ProjError {
    /// A conversion error
    #[error("The conversion failed with the following error: {0}")]
    Conversion(String),
    /// An error that occurs when a path string originating in PROJ can't be converted to a CString
    #[error("Couldn't create a raw pointer from the string")]
    Creation(#[from] std::ffi::NulError),
    #[error("Couldn't convert bytes from PROJ to UTF-8")]
    Utf8Error(#[from] std::str::Utf8Error),
    #[error("Couldn't convert number to f64")]
    FloatConversion,
}

#[derive(Error, Debug)]
pub enum ProjCreateError {
    #[error("A nul byte was found in the PROJ string definition or CRS argument: {0}")]
    ArgumentNulError(NulError),
    #[error("The underlying PROJ call failed: {0}")]
    ProjError(String),
}

pub(crate) struct Errno(pub libc::c_int);

impl Errno {
    /// Return the error message associated with the error number.
    pub fn message(&self, context: *mut PJ_CONTEXT) -> String {
        let ptr = unsafe { proj_sys::proj_context_errno_string(context, self.0) };
        if ptr.is_null() {
            panic!("PROJ did not supply an error")
        } else {
            unsafe { _string(ptr).expect("PROJ provided an invalid error string") }
        }
    }
}

/// Easily get a String from the external library
pub(crate) unsafe fn _string(raw_ptr: *const c_char) -> Result<String, str::Utf8Error> {
    assert!(!raw_ptr.is_null());
    let c_str = CStr::from_ptr(raw_ptr);
    Ok(str::from_utf8(c_str.to_bytes())?.to_string())
}

/// Look up an error message using the error code
fn error_message(code: c_int) -> Result<String, str::Utf8Error> {
    unsafe {
        let rv = proj_errno_string(code);
        _string(rv)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_convert() {
        let crs_from = "EPSG:7415";
        // Because we have a boundingVolume.box. For a boundingVolume.region we need 4979.
        let crs_to = "EPSG:4978";
        let transformer = Proj::new_known_crs(&crs_from, &crs_to, None).unwrap();
        let result = transformer.convert((85285.279, 446606.813, 10.0)).unwrap();
        println!("{:?}", result);
        // assert_relative_eq!(result.x() as f64, 3923215.044, epsilon = 1e-2);
        // assert_relative_eq!(result.y() as f64, 299940.760, epsilon = 1e-2);
        // assert_relative_eq!(result.z() as f64, 5003047.651, epsilon = 1e-2);

        // [0] = {f64} 85185.2799868164
        // [1] = {f64} 446506.81397216802
        // [2] = {f64} -15.333460330963135
        // [3] = {f64} 85385.2799868164
        // [4] = {f64} 446706.81397216802
        // [5] = {f64} 62.881539669036869

        // [0] = {f64} 3923286.6789069851
        // [1] = {f64} 299847.35238179943
        // [2] = {f64} 5002965.2671014387
        // [3] = {f64} 3923160.3138962006
        // [4] = {f64} 300035.46416343335
        // [5] = {f64} 5003151.7442537257
    }
}
