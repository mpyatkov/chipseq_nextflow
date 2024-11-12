This script simplifies the downloading of PDF files from the UCSC browser.

By default, it only requires the session name as input.

At the start, it scans the latest RESULTS* directory for *top25*.xlsx files
and includes them as an additional parameter in the qsub script.

The output consists of combined PDF files, where each coordinate region
from the *top25*.xlsx file corresponds to a separate PDF page. In addition
to the PDF files, the script also generates log files for debugging.

USAGE:

./run.sh session1 [session2] [session3] ...

The script supports multiple sessions, separated by spaces. If a session
name contains spaces, enclose it in quotes, for example: "some session name".
