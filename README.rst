qp-ivar-trim Qiita Plugin
==============================

|Build Status| |Coverage Status|

`Qiita <https://github.com/biocore/qiita/>`__ (canonically pronounced *cheetah*)
is an analysis environment for microbiome (and other "comparative -omics")
datasets.

This package includes the qp-fastp Qiita plugin.

How to test this package?
-------------------------
In order to test the qp-ivar-trim package, a local
installation of Qiita should be running in test mode on the address
`https://localhost:21174`, with the default test database created in Qiita's
test suite. Also, if Qiita is running with the default server SSL certificate,
you need to export the variable `QIITA_SERVER_CERT` in your environment, so the
Qiita Client can perform secure connections against the Qiita server:

.. code-block:: bash

    $ export QIITA_SERVER_CERT=<QIITA_INSTALL_PATH>/qiita_core/support_files/server.crt

Credits
-------

This plugin was created with `Cookiecutter <https://github.com/audreyr/cookiecutter>`__
and the `qiita-spots/qp-template-cookiecutter <https://github.com/qiita-spots/qp-template-cookiecutter>`__
project template.

.. |Build Status| image:: https://travis-ci.org/qiita-spots/qp-ivar-trim.png?branch=master
   :target: https://travis-ci.org/qiita-spots/qp-ivar-trim
.. |Coverage Status| image:: https://coveralls.io/repos/qiita-spots/qp-ivar-trim/badge.png?branch=master
   :target: https://coveralls.io/r/qiita-spots/qp-ivar-trim