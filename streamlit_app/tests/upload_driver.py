"""Drive the app's upload path headlessly: hand the page a chosen file, render it.

Since issue #64 the uploader is in the sidebar, and the data page loads whatever that
chooser has handed over through ``PENDING_UPLOAD_KEY`` rather than calling ``st.file_uploader``
itself. So this driver seeds that key instead of stubbing the widget — the same seam the
sidebar's own callback writes, which is what keeps this a check of the load path rather than
of a stub.
"""

import os
import sys
from pathlib import Path

APP = Path(os.environ["MAFIGATE_APP_DIR"])
sys.path.insert(0, str(APP))

import streamlit as st

MAF = Path(os.environ["MAFIGATE_UPLOAD_FILE"])


class UploadStub:
    """What st.file_uploader hands back: a name plus bytes via getvalue()."""

    name = MAF.name
    type = "maf"

    def getvalue(self):
        return MAF.read_bytes()


from components.sidebar import PENDING_UPLOAD_KEY  # noqa: E402
from config.pipeline_params import pipeline_params  # noqa: E402
from page_modules.data_loading import show_data_loading_page  # noqa: E402

if "filter_params" not in st.session_state:
    st.session_state.filter_params = pipeline_params("somatic")

# Re-seeded on every run, because AppTest instances in this process share one Streamlit
# session: the load stamps its token, and a second driver run pointed at a *different* MAF
# has to be seen as a different file rather than as the same one arriving again.
st.session_state[PENDING_UPLOAD_KEY] = UploadStub()

show_data_loading_page()

data = (st.session_state["maf_data"] if "maf_data" in st.session_state else None)
filtered = (st.session_state["filtered_data"] if "filtered_data" in st.session_state else None)
st.text(f"UPLOAD_ROWS={0 if data is None else len(data)}")
st.text(f"UPLOAD_PASSED={'none' if filtered is None else len(filtered)}")
st.text(f"UPLOAD_DP_DTYPE={'none' if data is None else data['DP'].dtype}")
