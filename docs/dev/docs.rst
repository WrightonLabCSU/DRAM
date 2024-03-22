Editing These Docs
=========================


Some projects keep their documentation in a separate wiki and some use text files and MarkDown files. These are good options but both have draw backs that have driven us to to pick ReadTheDocs for this project.

Read The Docs / Sphinx
------------------

Working with ReadTheDocs in many cases requires a few extra steps. There are two things you must know to understand why. First ReadTheDocs is providing a free service to our open-source projects, and only open source projects. We don’t always open our source code immediately and even when we do we often want to make updates that documentation and we don’t want to have to publish the docs mid edit to see our works. The disparity between the public read the docs service and the internal editing process requires that the docs be viewed locally.

The second thing that you need to know is what ReadTheDocs is doing for us on their servers so you can do the same locally. In short we write our docs in a markup language  called sphinx this is a lot like the mark down that we use on git hub but more powerful a term here meaning it can do more things like include code directly from the repository in the docs or include docs embeded in the code. ReadTheDocs provides the computation to compile  the sphinx code and host the html that is produced. But if you want to look at non public sphics docks you need to compile the sphix code so here is an example of how to do that.

Example:
-------

First download the DRAM2 repository, and if necessary checkout a branch containing work in progress.

.. code-block:: bash

   git clone https://github.com/WrightonLabCSU/DRAM2.git # You will need to enter a password
   cd DRAM2  # Chang to that dir

If you want to read docs before they are published you may be looking at a specific branch. If you are not looking into a branch, you can skip this step and work off of main.

.. code-block:: bash
   git checkout upgrade_docs_branch


You need to install the requirements with pip. You probably have pip installed but if you don't you can install it into a conda environment with 'conda install pip'. Note that if you install these in a conda environment, you need to use the same conda environment every time you compile the docs.

.. code-block:: bash

   conda install pip
   pip install -r docs/requirements.txt

Now that you are in the DRAM2 repo you just need to build the documentation. A good way to do that is with sphinx-autobuild.

You may need to install this tool, use pip to do that::

   pip install sphinx-autobuild


Then you can build::

    sphinx-autobuild docs docs/_build/html

You will see a message like `Serving on http://127.0.0.1:8000` use the link and you will be able to view the docs in all their glory.

[You can learn more about sphinx here](https://www.sphinx-doc.org/en/master/usage/index.html).

