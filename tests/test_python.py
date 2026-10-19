"""Tests for the kinegen Python bindings (run: python tests/test_python.py)."""
import math
import os
import sys

import numpy as np

import kinegen

HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES = os.path.join(HERE, "..", "examples")


def test_streaming_and_batches():
    prior = kinegen.load_prior(os.path.join(EXAMPLES, "randomized_tasks.yaml"))
    task = kinegen.sample_task(prior, seed=1, index=3)
    assert task.task_id.startswith("task_3_")
    assert "X.mass" in task.parameters
    det = kinegen.DetectorModel(task)
    truth = kinegen.generate_truth(task, 0)
    reco = det.apply(truth)
    assert reco.feature_names == task.feature_names
    for o in reco.objects:
        assert len(o.features) == len(task.feature_names)
        if o.truth_index >= 0:
            assert truth.particles[o.truth_index].object_type != "invisible"
    b = kinegen.generate_batch(task, 2000, start=0, threads=4)
    n, m, _ = b["p4"].shape
    assert n == 2000 and m == b["n_objects"].max()
    assert b["features"].shape == (n, m, len(task.feature_names))
    assert (b["mask"].sum(axis=1) == b["n_objects"]).all()
    assert (b["type"][~b["mask"]] == -1).all()
    assert np.isnan(b["features"][~b["mask"]]).all()
    # Batches are deterministic and independent of threading and batch boundaries.
    b1 = kinegen.generate_batch(task, 100, start=1000, threads=1)
    b2 = kinegen.generate_batch(task, 2000, start=0, threads=3)
    k = b1["n_objects"].max()
    assert np.array_equal(b1["p4"], b2["p4"][1000:1100, :k])
    assert b1["event_index"][0] == 1000
    # Truncation keeps the leading objects; n_objects keeps the full count.
    t = kinegen.generate_batch(task, 50, max_objects=2)
    assert t["p4"].shape[1] == 2


def test_environments_and_serialization():
    prior = kinegen.load_prior(os.path.join(EXAMPLES, "detector_randomization.yaml"))
    task = kinegen.sample_task(prior, seed=2)
    envs = [kinegen.sample_environment(prior, 5, i) for i in range(3)]
    assert len({e.parameters["detector.jet.resolution.stochastic"] for e in envs}) == 3
    truth = kinegen.generate_truth(task, 4)
    det = kinegen.DetectorModel(task)
    mets = [det.apply(truth, e).met for e in envs]
    assert len(set(mets)) == 3
    again = kinegen.TaskConfig.from_yaml(task.to_yaml())
    assert again.to_yaml() == task.to_yaml()
    assert kinegen.generate_truth(again, 4).to_yaml() == truth.to_yaml()
    env = kinegen.Environment.from_yaml(envs[1].to_yaml())
    assert det.apply(truth, env).to_yaml() == det.apply(truth, envs[1]).to_yaml()
    assert kinegen.assign_split(task, by="environment", env=envs[0]) in ("train", "val", "test")


def test_truth_ancestry_and_conservation():
    prior = kinegen.load_prior(os.path.join(EXAMPLES, "fixed_process.yaml"))
    task = kinegen.sample_task(prior)
    ev = kinegen.generate_truth(task, 0)
    incoming = np.sum([p.p4 for p in ev.particles if p.status == "incoming"], axis=0)
    final = np.sum([p.p4 for p in ev.particles if p.status in ("final", "recoil")], axis=0)
    assert np.allclose(incoming, final, rtol=1e-9, atol=1e-6)
    for i in ev.final_state():
        anc = ev.ancestors(i)
        if ev.particles[i].name != "recoil":
            assert ev.particles[anc[-1]].name in ("t", "tbar")


def test_invalid_prior_raises():
    try:
        kinegen.parse_prior("particles: {A: {mass: 10, object: jet, decays: [[A, A]]}}\nproduction: [A]\n")
    except ValueError as e:
        assert "cyclic" in str(e)
    else:
        raise AssertionError("expected ValueError")
    bad = "particles:\n  X: {pdg: 1, mass: 10, object: jet, decays: [[Y, Y]]}\n  Y: {pdg: 1, mass: 8, object: jet}\nproduction: [X]\n"
    try:
        kinegen.sample_task(kinegen.parse_prior(bad + "task: {validation: error}\n"))
    except kinegen.InvalidTask as e:
        assert "kinematically forbidden" in str(e)
    else:
        raise AssertionError("expected InvalidTask")


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print("[ PASS ]", name)
            except Exception as exc:  # noqa: BLE001
                failed += 1
                print("[ FAIL ]", name, repr(exc))
    sys.exit(1 if failed else 0)
