/* crate use */
use rand::SeedableRng;

/* project use */
use rustyread;

/* criterion use */
use criterion::criterion_group;

fn add_error_worker(
    target: f64,
    seq: &[u8],
    model: &rustyread::model::Error,
) -> rustyread::simulate::error::Changes {
    let mut changes = rustyread::simulate::error::Changes::new();
    let mut local_rng = rand::rngs::SmallRng::from_entropy();
    rustyread::simulate::error::add_error(
        7,
        target as f64,
        &seq,
        &mut changes,
        &model,
        &mut local_rng,
    );

    changes
}

fn add_error(c: &mut criterion::Criterion) {
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);

    let seq = rustyread::random_seq(100_000, &mut rng);
    let model = rustyread::model::Error::random(7);

    let mut g = c.benchmark_group("Add error");
    g.sampling_mode(criterion::SamplingMode::Flat);

    for target in (1000..=20_000).step_by(1000) {
        g.bench_with_input(
            criterion::BenchmarkId::new("target edit", target),
            &target,
            |b, &target| {
                b.iter(|| criterion::black_box(add_error_worker(target as f64, &seq, &model)))
            },
        );
    }
}

fn setup(c: &mut criterion::Criterion) {
    add_error(c);
}

criterion::criterion_group!(benches, setup);

criterion::criterion_main!(benches);
