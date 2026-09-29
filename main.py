from radmetrics.pipeline import run_pipeline
from radmetrics.io import run_bash_in_python

from pathlib import Path

base_dir = Path(__file__).resolve().parent
script_path = base_dir / "scripts" / "gathering.sh"

if __name__ == "__main__":
    run_bash_in_python(str(script_path))
    run_pipeline()

    print("=" * 22)
    print("* Successfully Done *")
    print("=" * 22)
