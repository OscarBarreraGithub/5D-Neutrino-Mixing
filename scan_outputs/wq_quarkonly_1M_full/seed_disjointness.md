# WQ Quark-Only 1M Seed Disjointness

- tasks: 50
- M_KK tiles/task: 10
- draws/M_KK/task: 2000
- tile seed stride: 1000003
- shard seed block: 20000000
- per-shard base: 202606040000 + 20000000 * task_id
- tile seed: base + 1000003 * tile_id, for tile_id=0..9
- draw seed: tile_seed + draw_idx
- required block > stride*n_tiles + n_draws = 10002030; chosen block = 20000000
- conclusion: draw-seed intervals are disjoint within and across all shards.
