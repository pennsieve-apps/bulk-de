"""
One image, two runtimes. Fixes the thread count from THREADS (or the cgroup CPU
quota, see safe_inputs.cpu_threads) before numpy is imported; then, as a Lambda
function (AWS_LAMBDA_RUNTIME_API is set), hands over to the Lambda runtime
interface client with processor.handler, or, as a container, runs main() once.
"""
import os
import sys

from processor.safe_inputs import cpu_threads

raw = os.environ.get("THREADS", "").strip()
n = cpu_threads(int(raw) if raw.isdigit() else 0)
for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[var] = str(n)
os.environ["BULK_DE_THREADS"] = str(n)
os.environ["BULK_DE_THREADS_SOURCE"] = "THREADS" if raw.isdigit() and int(raw) > 0 else "CPU quota"

if os.environ.get("AWS_LAMBDA_RUNTIME_API"):
    os.execvp(sys.executable, [sys.executable, "-m", "awslambdaric", "processor.handler.handler"])

from processor.main import main  # noqa: E402

sys.exit(main())
