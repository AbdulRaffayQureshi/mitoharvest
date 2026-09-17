import time
import logging
import http.client
import socket
from urllib.error import HTTPError, URLError
from Bio import Entrez
import config

RETRYABLE_ERRORS = (
    HTTPError, URLError, RuntimeError,
    http.client.HTTPException, ConnectionError, socket.timeout, TimeoutError,
)

log = logging.getLogger("ncbi_client")

Entrez.email = config.NCBI_EMAIL
if config.NCBI_API_KEY:
    Entrez.api_key = config.NCBI_API_KEY


class NCBIClient:
    def __init__(self, delay=config.REQUEST_DELAY):
        self.delay = delay
        self._last_call = 0.0

    def _throttle(self):
        elapsed = time.time() - self._last_call
        wait = self.delay - elapsed
        if wait > 0:
            time.sleep(wait)
        self._last_call = time.time()

    def _call(self, func, parser=lambda h: h.read(), **kwargs):
        last_error = None
        for attempt in range(1, config.MAX_RETRIES + 1):
            self._throttle()
            try:
                handle = func(**kwargs)
                result = parser(handle)
                handle.close()
                return result
            except RETRYABLE_ERRORS as e:
                last_error = e
                wait = config.RETRY_BACKOFF_BASE ** attempt
                log.warning(f"NCBI call failed ({e}), retry {attempt}/{config.MAX_RETRIES} in {wait:.1f}s")
                time.sleep(wait)
        raise RuntimeError(f"NCBI call failed after {config.MAX_RETRIES} retries: {last_error}")

    def esearch(self, query, retstart=0, retmax=config.SEARCH_PAGE_SIZE, db="nucleotide"):
        return self._call(
            Entrez.esearch, parser=Entrez.read,
            db=db, term=query, retstart=retstart, retmax=retmax,
        )

    def efetch_text(self, ids, rettype, retmode="text", db="nucleotide"):
        ids_str = ",".join(ids)
        return self._call(
            Entrez.efetch, parser=lambda h: h.read(),
            db=db, id=ids_str, rettype=rettype, retmode=retmode,
        )

    def esummary(self, ids, db="pubmed"):
        ids_str = ",".join(str(i) for i in ids)
        return self._call(
            Entrez.esummary, parser=Entrez.read,
            db=db, id=ids_str,
        )
# import time
# import logging
# import http.client
# import socket
# from urllib.error import HTTPError, URLError
# from Bio import Entrez
# import config

# RETRYABLE_ERRORS = (
#     HTTPError, URLError, RuntimeError,
#     http.client.HTTPException, ConnectionError, socket.timeout, TimeoutError,
# )

# log = logging.getLogger("ncbi_client")

# Entrez.email = config.NCBI_EMAIL
# if config.NCBI_API_KEY:
#     Entrez.api_key = config.NCBI_API_KEY


# class NCBIClient:
#     def __init__(self, delay=config.REQUEST_DELAY):
#         self.delay = delay
#         self._last_call = 0.0

#     def _throttle(self):
#         elapsed = time.time() - self._last_call
#         wait = self.delay - elapsed
#         if wait > 0:
#             time.sleep(wait)
#         self._last_call = time.time()

#     def _call(self, func, parser=lambda h: h.read(), **kwargs):
#         last_error = None
#         for attempt in range(1, config.MAX_RETRIES + 1):
#             self._throttle()
#             try:
#                 handle = func(**kwargs)
#                 result = parser(handle)
#                 handle.close()
#                 return result
#             except RETRYABLE_ERRORS as e:
#                 last_error = e
#                 wait = config.RETRY_BACKOFF_BASE ** attempt
#                 log.warning(f"NCBI call failed ({e}), retry {attempt}/{config.MAX_RETRIES} in {wait:.1f}s")
#                 time.sleep(wait)
#         raise RuntimeError(f"NCBI call failed after {config.MAX_RETRIES} retries: {last_error}")

#     def esearch(self, query, retstart=0, retmax=config.SEARCH_PAGE_SIZE):
#         return self._call(
#             Entrez.esearch, parser=Entrez.read,
#             db="nucleotide", term=query, retstart=retstart, retmax=retmax,
#         )

#     def efetch_text(self, ids, rettype, retmode="text"):
#         ids_str = ",".join(ids)
#         return self._call(
#             Entrez.efetch, parser=lambda h: h.read(),
#             db="nucleotide", id=ids_str, rettype=rettype, retmode=retmode,
#         )