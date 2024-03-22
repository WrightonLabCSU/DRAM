.. RMNP_Pipeline documentation master file, created by
   sphinx-quickstart on Thu Jul 14 14:52:41 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introducing Tree Kit For DRAM
=========================================


.. toctree::
   :maxdepth: 2
   :caption: Contents:

## Usage
Coming soon

### Input Options

### Output to expect

## Status
At this point there are 2 tree files (nxr_nar, and pmoa_amoa) and one working tree (nxr_nar). Each has a mapping that takes leaves and reduces them to calls, then the call that is closessed is applied. but all data is saved a stronger tree based method is still in dev and may not be needed? These mappings are from the excell file you provied but modified it. To start with the calls did not make sense to me, nore did the adjectives, so I droped them. I droped the genome info keeping only the gene, the call+adjective and the tree. then I split them by tree. 

At this time there are 4 mutualy exclusive clasifications in the Nxr/Nar tree and they are 'other-None', 'nxr-None', 'nxr-Nitrifier', 'nxr/nar-N utilization', 'narG-N reducer', this is the first thing that needs to be tweeked


The Process
=====================

Below you can see a more complex overview of annotation process. 


.. raw:: html
   <div class="mxgraph" style="max-width:100%;border:1px solid transparent;" data-mxgraph="{&quot;highlight&quot;:&quot;#0000ff&quot;,&quot;nav&quot;:true,&quot;resize&quot;:true,&quot;toolbar&quot;:&quot;zoom layers tags lightbox&quot;,&quot;edit&quot;:&quot;_blank&quot;,&quot;url&quot;:&quot;https://raw.githubusercontent.com/rmFlynn/DRAM_Trees/main/docs/figs/tree_kit_fig.drawio&quot;}"></div>
   <script type="text/javascript" src="https://viewer.diagrams.net/embed2.js?&fetch=https%3A%2F%2Fraw.githubusercontent.com%2FrmFlynn%2FDRAM_Trees%2Fmain%2Fdocs%2Ffigs%2Ftree_kit_fig.drawio"></script>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
