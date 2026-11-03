import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from classify import classify_call, ribosomal_fail


def test_called():
    assert classify_call(1.2, True, 10) == "called_up"
    assert classify_call(-0.5, True, 99) == "called_down"

def test_boundaries():
    assert classify_call(0.1, False, 74.99) == "not_called_small"
    assert classify_call(0.1, False, 75) == "not_called_mid"
    assert classify_call(0.1, False, 94.99) == "not_called_mid"
    assert classify_call(-2.0, False, 95) == "not_called_large"

def test_ribosomal_median():
    assert ribosomal_fail([80, 80, 76], [1, -1, 1]) is True      # median 80 > 75
    assert ribosomal_fail([75, 75, 75], [1, -1, 1]) is False     # median 75 not > 75

def test_ribosomal_sign():
    # 5 of 5 positive = 100% > 80%
    assert ribosomal_fail([10] * 5, [1] * 5) is True
    # 4 of 5 positive = 80%, not > 80%
    assert ribosomal_fail([10] * 5, [1, 1, 1, 1, -1]) is False
    # 9 of 10 negative = 90%
    assert ribosomal_fail([10] * 10, [-1] * 9 + [1]) is True
