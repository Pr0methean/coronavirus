import filecmp
import os
import unittest

import redis
from funcy import chunks

from design_guides import CORONAVIRUS_CONSENSUS, HOST_LUNG_TISSUE, main, \
    open_as_mmap


class EndToEndTest(unittest.TestCase):

    def test_end_to_end(self):
        test_redis = redis.Redis(host='localhost', port=6379, db='15')
        test_redis.flushdb()
        reused_args = {
            'r': test_redis,
            'host_path': os.path.join("host", HOST_LUNG_TISSUE),
            'target_path': os.path.join("alignments", CORONAVIRUS_CONSENSUS),
            'vectors_path': os.path.join("tmp", "test_vectors"),
            'index_path': os.path.join("tmp", "test_index.pkl"),
            'out_path': os.path.join("tmp", "test_output"),
            'target_id': 'nCoV',
            'k': 28,
            'max_mismatches_for_side_effect': 8
        }
        try:
            main(rebuild_index=True,
                 use_existing_vectors=False,
                 **reused_args)
            self.verify_outputs()
            # Cover branches where index is rebuilt but existing vectors are
            # used
            main(rebuild_index=True,
                 use_existing_vectors=True,
                 **reused_args)
            self.verify_outputs()
            # Cover branches where existing index is used
            main(rebuild_index=False,
                 **reused_args)
            self.verify_outputs()
        finally:
            test_redis.flushdb()

    def verify_outputs(self):
        self.assertTrue(
            filecmp.cmp(os.path.join("testdata", "lung_tissue_output"),
                        os.path.join("tmp", "test_output"),
                        shallow=False),
            f"File test_output in tmp doesn't match reference file "
            f"lung_tissue_output in testdata")
        self.assertEqual(
            {b''.join(bytelist) for bytelist in
             chunks(28, 28, open_as_mmap(
                 os.path.join("testdata", "lung_tissue_vectors")))},
            {b''.join(bytelist) for bytelist in
             chunks(28, 28, open_as_mmap(os.path.join("tmp", "test_vectors")))}
        )
        # Doesn't work because BallTree doesn't implement __eq__:
        # self.assertEqual(load_from_pickle(os.path.join("testdata",
        # "lung_tissue_index.pkl")),
        #                  load_from_pickle(os.path.join("tmp",
        #                  "test_index.pkl")))


if __name__ == '__main__':
    unittest.main()
