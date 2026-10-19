"""Feed randomised tasks directly into ML training (no intermediate files).

    PYTHONPATH=build python examples/train_feed.py

Each training "episode" is one task: a support set and a query set of events
drawn from the same sampled physics task, observed through a detector
environment that is itself randomised. Tasks are split into train/val/test by a
deterministic hash of their topology, so held-out topologies never leak.
"""
import os

import numpy as np

import kinegen

HERE = os.path.dirname(os.path.abspath(__file__))


def episodes(prior, split="train", seed=0, n_support=256, n_query=256, threads=4):
    index = 0
    while True:
        task = kinegen.sample_task(prior, seed=seed, index=index)
        index += 1
        if kinegen.assign_split(task, by="topology") != split:
            continue
        # A fresh detector environment per episode (domain randomisation).
        env = kinegen.sample_environment(prior, seed + 1, index)
        support = kinegen.generate_batch(task, n_support, start=0, env=env, threads=threads)
        query = kinegen.generate_batch(task, n_query, start=n_support, env=env, threads=threads)
        yield task, env, support, query


def main():
    prior = kinegen.load_prior(os.path.join(HERE, "randomized_tasks.yaml"))
    for k, (task, env, support, query) in enumerate(episodes(prior)):
        x = support["ptetaphim"]                # (n, max_objects, 4)
        feats = support["features"]             # (n, max_objects, n_features), NaN = missing
        mask = support["mask"]                  # (n, max_objects) valid objects
        label = support["label"]                # (n,) mixture component of each event
        print(f"{task.task_id}: {task.topology}")
        print(f"  env {env.env_id}: jet stochastic term {env.parameters['detector.jet.resolution.stochastic']:.2f}")
        print(f"  objects/event {mask.sum(1).mean():.2f}, features {support['feature_names']}")
        print(f"  label balance {np.bincount(label)}, query shape {query['ptetaphim'].shape}")
        if k == 2:
            break


if __name__ == "__main__":
    main()
