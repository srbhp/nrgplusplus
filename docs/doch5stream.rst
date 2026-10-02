HDF5 I/O
========

The ``h5stream`` utilities provide a compact file-based interface for storing and loading HDF5 datasets
used in NRG workflows, including eigenvalues, operator data, and saved iteration results.

How the wrapper maps data to HDF5
--------------------------------

``h5stream::h5stream`` owns an HDF5 file handle and exposes typed read and write helpers. Its constructor opens
the file in one of four modes: ``r`` for read-only, ``rw`` for read/write, ``x`` for exclusive creation, or
``tr`` to create or truncate the file. The default mode is ``tr``, so opening an existing path with the
default will replace its contents.

Before writing, ``get_datatype_for_hdf5`` maps supported C++ scalar types to their native HDF5 types. A vector
is stored as a one-dimensional dataset; nested vectors are written as a sequence of datasets with numeric
suffixes appended to the base name. The matching nested-vector reader loads that sequence in order. Raw
pointer-and-size overloads support callers that already manage contiguous storage.

Metadata is stored as HDF5 attributes. The ``dspace`` and ``gspace`` helpers bind attribute operations to a
dataset or group, while ``getDataspace``, ``getGroup``, and ``createGroup`` provide the corresponding file
navigation operations. Call ``close`` when file access is complete.

For the generated class reference, see the ``h5stream::h5stream`` `class documentation <api/classh5stream_1_1h5stream.html>`_.

