{
  description = "FiniteFields: Julia, Lean 4, and Quartz development environment";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/55bfb4d7437bd4e09ab7cba58a43bd68eb310044";

  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" "aarch64-darwin" ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems
        (system: f nixpkgs.legacyPackages.${system});

      useSystemJulia = pkgs: ''
        if command -v julia >/dev/null 2>&1; then
          echo "FiniteFields: using system Julia ($(command -v julia))"
        else
          PATH="${pkgs.julia-bin}/bin:$PATH"
          echo "FiniteFields: no system Julia found, using Nix $(julia --version)"
        fi
      '';

      useSystemLean = pkgs: ''
        if command -v elan >/dev/null 2>&1 || command -v lake >/dev/null 2>&1; then
          echo "FiniteFields: using system Lean toolchain"
        else
          PATH="${pkgs.elan}/bin:$PATH"
          echo "FiniteFields: no system Lean found, using Nix elan"
        fi
      '';
    in
    {
      devShells = forAllSystems (pkgs: {
        julia = pkgs.mkShell {
          packages = [ pkgs.git pkgs.nodejs_22 ];
          shellHook = useSystemJulia pkgs;
        };

        lean = pkgs.mkShell {
          packages = [ pkgs.git ];
          shellHook = useSystemLean pkgs;
        };

        quartz = pkgs.mkShell {
          packages = [ pkgs.git pkgs.nodejs_22 ];
        };

        default = pkgs.mkShell {
          packages = [ pkgs.git pkgs.nodejs_22 ];
          shellHook = useSystemJulia pkgs + useSystemLean pkgs;
        };
      });
    };
}
