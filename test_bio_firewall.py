import csv
import os

from redis import Redis

from bio_firewall import count_records, get_kmers, make_hosts, make_targets, \
    predict_side_effects

TEST_REDIS_ARGS = {'db': 15}


def test_count_records():
    path = os.path.join("data", "test", "unit_test_host.fa")
    total = count_records(path)
    assert total == 2


def test_get_kmers():
    assert list(get_kmers('actg', k=2, stringify=0)) == ['ac', 'ct', 'tg']
    wc_kmers = list(get_kmers('wc', k=2, stringify=0))
    for kmer in ['ac', 'tc']:
        assert kmer in wc_kmers
    nt_kmers = list(get_kmers('nt', k=2, stringify=0))
    for kmer in ['at', 'ct', 'tt', 'gt']:
        assert kmer in nt_kmers
    wn_kmers = list(get_kmers("wn", k=2, stringify=0))
    for kmer in ['aa', 'at', 'ac', 'ag', 'ta', 'tt', 'tc', 'tg']:
        assert kmer in wn_kmers


def test_make_hosts():
    r, k = Redis(**TEST_REDIS_ARGS), 28
    r.flushdb()
    path = os.path.join("data", "test", "unit_test_host.fa")
    hosts = make_hosts(path=path, k=k, redis_args=TEST_REDIS_ARGS)
    hosts_snapshot_path = os.path.join("data", "snapshots", "hosts")
    with open(hosts_snapshot_path, "r") as hosts_snapshot_file:
        hosts_snapshot = hosts_snapshot_file.read().splitlines()
    print("hosts", hosts)
    print("hosts_snapshot", hosts_snapshot)
    for host in hosts:
        assert host in hosts_snapshot


def test_make_targets():
    r, k = Redis(**TEST_REDIS_ARGS), 28
    r.flushdb()
    test_alignment_path = os.path.join("data", "test", "unit_test_target.clu")
    targets = make_targets(path=test_alignment_path, id="nCoV", k=k, db=r)
    targets_snapshot_path = os.path.join("data", "snapshots", "targets")
    with open(targets_snapshot_path, "r") as targets_snapshot_file:
        targets_snapshot = list(csv.reader(targets_snapshot_file))
    # print("targets", targets)
    # print("targets_snapshot", targets_snapshot)
    for target, snapshot in zip(targets, targets_snapshot):
        assert target[0] == snapshot[0]
        assert str(target[1]).strip(' ') == snapshot[1].strip(' ')


def test_predict_side_effects():
    k, cutoff = 28, 1
    r = Redis(**TEST_REDIS_ARGS)
    r.flushdb()
    test_host_path = os.path.join("data", "test", "unit_test_host.fa")
    make_hosts(path=test_host_path, k=k, redis_args=TEST_REDIS_ARGS)
    test_alignment_path = os.path.join("data", "test", "unit_test_target.clu")
    make_targets(path=test_alignment_path, id="nCoV", k=k, db=r)
    good_targets_snapshot_path = os.path.join("data", "snapshots",
                                              "good_targets")
    with open(good_targets_snapshot_path, "r") as good_targets_snapshot_file:
        good_targets_snapshot = list(csv.reader(good_targets_snapshot_file))
    good_targets = predict_side_effects(k=k, cutoff=cutoff, db=r)
    print("good_targets", good_targets)
    print("good_targets_snapshot", good_targets_snapshot)
    for good_target, snapshot in zip(good_targets, good_targets_snapshot):
        assert good_target[0] == snapshot[0]
        assert str(good_target[1]).strip(' ') == snapshot[1].strip(' ')


# save a sample of successful data to compare with future results (inspect!)
def snapshot():
    cutoff = 1
    k = 28
    r = Redis(**TEST_REDIS_ARGS)
    r.flushdb()
    host_path = os.path.join("data", "test", "unit_test_host.fa")
    hosts_results_path = os.path.join("data", "snapshots", "hosts")
    hosts = make_hosts(path=host_path, k=k)
    print("hosts:")
    print(hosts)
    print("")
    with open(hosts_results_path, "w+") as hosts_results_file:
        hosts_results_file.writelines([f"{h}\n" for h in hosts])
    alignment_path = os.path.join("data", "test", "unit_test_target.clu")
    targets_results_path = os.path.join("data", "snapshots", "targets")
    targets = make_targets(path=alignment_path, k=k)
    print("targets:")
    print(targets)
    print("")
    with open(targets_results_path, "w+") as targets_results_file:
        targets_results_file.writelines([f"{t[0]}, {t[1]}\n" for t in targets])
    good_targets = predict_side_effects(k=k, cutoff=cutoff)
    good_targets_results_path = os.path.join("data", "snapshots",
                                             "good_targets")
    print("good_targets:")
    print(good_targets)
    print("")
    with open(good_targets_results_path, "w+") as good_targets_results_file:
        good_targets_results_file.writelines(
            [f"{gt[0]}, {gt[1]}\n" for gt in good_targets])

# snapshot()

