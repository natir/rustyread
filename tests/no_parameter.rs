mod common;

#[test]
#[should_panic]
fn no_parameter() {
    common::init();

    common::run_process("./target/debug/rustyread", &["-vvvvv"]);
}
