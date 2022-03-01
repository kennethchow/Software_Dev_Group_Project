from asyncio import ThreadedChildWatcher
from concurrent.futures import thread
from website import create_app
import os
# ======= App Creation (from __init__.py ======= #
app = create_app()

# ======= Run File ======= #
if __name__ == '__main__':
    app.run(debug=True,
            threaded=True,
            host = '0.0.0.0',
            port = int(os.environ.get("PORT", 8080)))

