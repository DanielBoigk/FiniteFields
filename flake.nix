{
  description = "FiniteFields: Julia, Lean 4, and Quartz development environment";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/55bfb4d7437bd4e09ab7cba58a43bd68eb310044";

  # Toolchains are fetched lazily. A shell's hook looks for the tool on the system first and
  # only builds (downloads) the Nix copy when it is missing. The hooks deliberately do not
  # interpolate `${pkgs.julia-bin}` or `${pkgs.elan}`: a store path in the hook string would
  # make the shell depend on it, and Nix would download it even when the system has the tool.
  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" "x86_64-darwin" "aarch64-darwin" ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems
        (system: f system nixpkgs.legacyPackages.${system});

      # `probe`: commands whose presence means the system already has the tool.
      # `attr`: the nixpkgs attribute to build when it does not. It is built from the pinned
      # nixpkgs source already in the store (`${nixpkgs}`), so nothing is re-locked or re-fetched
      # besides the tool itself.
      lazyTool = system: { label, probe, attr }: ''
        if ${nixpkgs.lib.concatMapStringsSep " || " (c: "command -v ${c} >/dev/null 2>&1") probe}; then
          echo "FiniteFields: using system ${label} ($(command -v ${builtins.head probe}))"
        else
          echo "FiniteFields: no system ${label} found, fetching it with Nix ..."
          if _ff_out="$(nix build --no-link --print-out-paths "path:${nixpkgs}#legacyPackages.${system}.${attr}")"; then
            PATH="$_ff_out/bin:$PATH"
            echo "FiniteFields: using Nix ${label} ($_ff_out)"
          else
            echo "FiniteFields: could not fetch ${label} with Nix" >&2
          fi
          unset _ff_out
        fi
      '';

      julia = system: lazyTool system { label = "Julia"; probe = [ "julia" ]; attr = "julia-bin"; };
      lean = system: lazyTool system { label = "Lean toolchain (elan)"; probe = [ "elan" "lake" ]; attr = "elan"; };
    in
    {
      devShells = forAllSystems (system: pkgs: {
        # Julia work. No Lean tooling, no Node.
        julia = pkgs.mkShell {
          packages = [ pkgs.git ];
          shellHook = julia system;
        };

        # Lean work. elan is only ever fetched from this shell.
        lean = pkgs.mkShell {
          packages = [ pkgs.git ];
          shellHook = lean system;
        };

        # Vault and site work: Node.js 22 for Quartz.
        quartz = pkgs.mkShell {
          packages = [ pkgs.git pkgs.nodejs_22 ];
        };

        # Default: Julia plus the Quartz tools. Lean stays opt-in (`nix develop .#lean`).
        default = pkgs.mkShell {
          packages = [ pkgs.git pkgs.nodejs_22 ];
          shellHook = julia system;
        };
      });
    };
}
