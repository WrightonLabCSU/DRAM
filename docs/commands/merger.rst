Merge data
==========

Many workflows benefit from running DRAM2 processes separately in a map-reduce or scatter and gather method. This is especially the case with Snakemake. When this is the case you then need to combine the results after the run completes. This is where the merge command comes in.

Limitations
^^^^^^^^^^^

At this time, only the merging of gene calls and annotations are supported. These are the first two steps and are the most time consuming parts of DRAM. It it almost always the case that once these two steps are run separately and merged creating a new distillate and annotation file is not a difficult task. However, as data sets get bigger and bigger DRAM2's post processing steps become more and more cumbersome. It is the long term goal that all DRAM2 outputs will be merged in the future.


.. raw:: html

   <body><div class="mxgraph" style="max-width:100%;border:1px solid transparent;" data-mxgraph="{&quot;highlight&quot;:&quot;#0000ff&quot;,&quot;nav&quot;:true,&quot;resize&quot;:true,&quot;toolbar&quot;:&quot;zoom layers tags lightbox&quot;,&quot;edit&quot;:&quot;_blank&quot;,&quot;xml&quot;:&quot;&lt;mxfile host=\&quot;app.diagrams.net\&quot; modified=\&quot;2023-03-13T20:58:46.082Z\&quot; agent=\&quot;5.0 (X11)\&quot; etag=\&quot;jnAd7CkboL1H8Lhm2q-L\&quot; version=\&quot;21.0.6\&quot; type=\&quot;github\&quot;&gt;&lt;diagram name=\&quot;Page-1\&quot; id=\&quot;3VVjTzhk1s7mqIIuyEGd\&quot;&gt;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&lt;/diagram&gt;&lt;/mxfile&gt;&quot;}"></div>
   <script type="text/javascript" src="https://viewer.diagrams.net/js/viewer-static.min.js"></script>

Inputs:
------

.. list-table::
   :widths: 25 25 50
   :header-rows: 1

   * - Input (Use)
     - Type
     - Notes
   * - DRAM2 Output Folders(ARGUMENT)
     - Directory Path, for DRAM Outputs
     - Unless you pass the force flag, the folders must contain DRAM2 project config files. If those conditions are not met, DRAM2 will look for the "genes" folder and the "raw.tsv" file by name and skip folders where neither are found.
   * - annotations file (-a/--annotations)
     - File path, pointing to an raw.tsv file
     - Use with the force flag. If you have your annotations in a custom location then you can combine them with this option. You can use it as many times as you need to point to as many files as you want.
   * - genes folder(-g/--genes)
     - Directory Path, for DRAM2 called genes
     - Use with the force flag. If you moved the genes directories for any reason, you can point to each of them by using this flag multipal times. The genes will be combind in the output if there are no colisions in the nameing. 
   * - force (-f/--force)
     - Flag
     - Skip any config checks and merge the files no mater the cost. If you use this command you except responsibility for the state of the result. Overwriten data, incompatible annotations may result. 


Output:
-------

The output will be a DRAM2 output directory as specified by the name passed to the :ref: `dram2` command's `-o/--output` option. This output directory will have all the data specified from the input options above and a new project_metadata file.  

examples:
----------

