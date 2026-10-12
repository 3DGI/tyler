{
  description = "Nix build for tyler";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";

  outputs =
    { self, nixpkgs, ... }:
    let
      supportedSystems = [
        "aarch64-darwin"
        "aarch64-linux"
        "x86_64-darwin"
        "x86_64-linux"
      ];

      forAllSystems = nixpkgs.lib.genAttrs supportedSystems;
    in
    {
      packages = forAllSystems (
        system:
        let
          pkgs = import nixpkgs { inherit system; };
          projData = pkgs.runCommand "tyler-proj-data" { } ''
            mkdir -p "$out/share/proj"
            cp -a ${pkgs.proj}/share/proj/. "$out/share/proj/"
            chmod -R u+w "$out/share/proj"
            if [ ! -e "$out/share/proj/nl_nsgi_nlgeo2018.tif" ]; then
              cp ${pkgs.fetchurl {
                url = "https://cdn.proj.org/nl_nsgi_nlgeo2018.tif";
                hash = "sha256-+OMsVr+JQPw/77wOQT60VGYz7VmLZpxjhWzd+DKJksA=";
              }} "$out/share/proj/nl_nsgi_nlgeo2018.tif"
            fi
            if [ ! -e "$out/share/proj/nl_nsgi_rdtrans2018.tif" ]; then
              cp ${pkgs.fetchurl {
                url = "https://cdn.proj.org/nl_nsgi_rdtrans2018.tif";
                hash = "sha256-dlODEZG0JOcVqQZGiWL8YAcc+3GjGGsqWPCYuri/Qd4=";
              }} "$out/share/proj/nl_nsgi_rdtrans2018.tif"
            fi
          '';
          nativeBuildInputs = [
            pkgs.cmake
            pkgs.llvmPackages.clang
            pkgs.makeWrapper
            pkgs.pkg-config
          ];
          buildInputs = [
            pkgs.libtiff
            pkgs.proj
            pkgs.sqlite
          ];
        in
        rec {
          default = tyler;

          tyler = pkgs.rustPlatform.buildRustPackage {
            pname = "tyler";
            version = (builtins.fromTOML (builtins.readFile ./Cargo.toml)).package.version;

            src = self;
            cargoLock = {
              lockFile = ./Cargo.lock;
              outputHashes = {
                "cityjson-index-0.10.0" = "sha256-1jq/5H4faqIqVJml6jM4+KrCuyiaxbSAIkEzY2krffg=";
                "meshopt-0.7.0" = "sha256-q9K1ZMKPw4RlCFoHpQ0ZrWrQ/h3WeP8TIcmEF8mopo4=";
              };
            };

            inherit nativeBuildInputs buildInputs;

            cargoBuildFlags = [ "--workspace" ];

            LIBCLANG_PATH = "${pkgs.llvmPackages.libclang.lib}/lib";
            PROJ_DATA = "${projData}/share/proj";

            doCheck = false;

            postInstall = ''
              wrapProgram "$out/bin/tyler" \
                --set-default PROJ_DATA "${projData}/share/proj"
              wrapProgram "$out/bin/cjconvert" \
                --set-default PROJ_DATA "${projData}/share/proj"
            '';
          };
        }
      );

      devShells = forAllSystems (
        system:
        let
          pkgs = import nixpkgs { inherit system; };
          projData = pkgs.runCommand "tyler-proj-data" { } ''
            mkdir -p "$out/share/proj"
            cp -a ${pkgs.proj}/share/proj/. "$out/share/proj/"
            chmod -R u+w "$out/share/proj"
            if [ ! -e "$out/share/proj/nl_nsgi_nlgeo2018.tif" ]; then
              cp ${pkgs.fetchurl {
                url = "https://cdn.proj.org/nl_nsgi_nlgeo2018.tif";
                hash = "sha256-+OMsVr+JQPw/77wOQT60VGYz7VmLZpxjhWzd+DKJksA=";
              }} "$out/share/proj/nl_nsgi_nlgeo2018.tif"
            fi
            if [ ! -e "$out/share/proj/nl_nsgi_rdtrans2018.tif" ]; then
              cp ${pkgs.fetchurl {
                url = "https://cdn.proj.org/nl_nsgi_rdtrans2018.tif";
                hash = "sha256-dlODEZG0JOcVqQZGiWL8YAcc+3GjGGsqWPCYuri/Qd4=";
              }} "$out/share/proj/nl_nsgi_rdtrans2018.tif"
            fi
          '';
        in
        {
          default = pkgs.mkShell {
            packages = [
              pkgs.cargo
              pkgs.clippy
              pkgs.cmake
              pkgs.llvmPackages.clang
              pkgs.pkg-config
              pkgs.rustc
              pkgs.rustfmt
	      pkgs.perf
	      pkgs.valgrind
            ];

            buildInputs = [
              pkgs.libtiff
              pkgs.proj
              pkgs.sqlite
            ];

            LIBCLANG_PATH = "${pkgs.llvmPackages.libclang.lib}/lib";
            PROJ_DATA = "${projData}/share/proj";
          };
        }
      );
    };
}
