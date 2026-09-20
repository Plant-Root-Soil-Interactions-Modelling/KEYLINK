"""
KEYLINK test package.

All model modules read their parameter files relative to the working directory, so this
package makes the model directory (the parent of tests/) the working directory and puts it
on sys.path before any test module imports keylink_core / keylink_functions. A non-interactive
matplotlib backend is selected because keylink_core imports pyplot.

Run everything with   python run_tests.py   (or: python -m unittest discover -s tests -t .,
or: pytest tests) from the model directory.
"""
import os
import sys

MODEL_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
os.environ.setdefault('MPLBACKEND', 'Agg')
if MODEL_DIR not in sys.path:
    sys.path.insert(0, MODEL_DIR)
os.chdir(MODEL_DIR)
