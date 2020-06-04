with import <nixpkgs> {};

mkShell {
  nativeBuildInputs = [
    cargo
    rustc
    rustfmt
  ];

  shellHook = ''
    export RUST_BACKTRACE=full
  '';
}
