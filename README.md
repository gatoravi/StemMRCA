# StemMRCA

Add **missing taxa to a phylogenetic tree** using taxonomic constraints, and
generate many random replicate trees (2011, Burleigh Lab, University of
Florida).

## Background

Large phylogenies built from sequence data almost always leave out species
that have no usable sequences. Comparative analyses such as diversification
rates and trait evolution can be biased if those species are simply dropped.
A common workaround is to add each missing species at a random position inside
the clade that taxonomy says it belongs to. Repeating this many times gives a
set of trees that captures the uncertainty in where the missing species
actually sit.

## What it does

For every missing taxon, you name two species already in the tree that bound
its group (for example, two members of its genus). StemMRCA:

1. finds the **most recent common ancestor (MRCA)** of those two species,
2. collects every branch in the clade below that MRCA, **including the stem
   branch** leading to it,
3. picks a branch at random with probability proportional to branch length,
   and attaches the new taxon there.

The process is repeated to produce the requested number of replicate trees in Newick format.

## Build

Requires a C++ compiler and the Boost headers (bundled under `include/`).

```bash
make        # produces ./StemAdd
```

## Usage

```bash
./StemAdd initial_tree.nwk taxa_to_add.txt 100 output_trees.nwk
```

| Argument | Meaning |
| --- | --- |
| `initial_tree.nwk` | Starting tree in Newick format, with branch lengths |
| `taxa_to_add.txt` | One line per missing taxon: `new_taxon  anchor_taxon_1  anchor_taxon_2` |
| `100` | Number of replicate trees to generate |
| `output_trees.nwk` | Output file of replicate trees |

Example inputs, including an antbird (Thamnophilidae) dataset and simulated
trees, are in [`datasets2/`](datasets2).

## Credits

Built on the tree library by André Wehe, Mukul Bansal and Oliver Eulenstein
(GPL-2.0), used for tree I/O and LCA queries.
