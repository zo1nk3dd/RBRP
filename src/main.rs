use rbrp::model::*;

fn main() {
    let mut model = CompactModel::new("instances/03-Bari10_RBRP.txt", 1, Some(0.1));
    model.solve();
}
