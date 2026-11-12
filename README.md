# Common-Neighbor-Count-Based Representative Possible World Finding on Uncertain Graphs

## 1 Introduction

This is a description of the code used for the experiments in the paper entitled *Common-Neighbor-Count-Based Representative Possible World Finding on Uncertain Graphs*.

Given an uncertain graph with probabilistic edges, a representative possible world (RPW) is a deterministic graph derived from the uncertain graph where a designated structural feature closely approximates its expected value. We study the Common-neighbor-count-based Representative Possible World (CRPW) problem, which seeks the possible world that best preserves the expected numbers of common neighbors between node pairs.

We provide three methods, as shown in Table 1.

**Table 1. Our Methods**

| ALGORITHM | DESCRIPTION |
| --- | --- |
| P+RSR | The basic algorithm, consisting of probability-based initialization (P), which produces an initial possible world, and random selection-based refinement (RSR), which reduces its total discrepancy. |
| P+CR | P+CR optimizes the efficiency of RSR with integer-based refinement (CR), which replaces the floating-point computation in each iteration with integer counting. |
| P+BCR | P+BCR builds upon CR with Beta-based integer refinement (BCR), which adaptively terminates the refinement instead of running a fixed number of iterations. |

The runner also includes the ADR, TRPW, and ABM baselines. For inspection, it saves the possible world produced by probability-based initialization as PIN. In Python identifiers and command-line values, the `+` in the paper names is written as an underscore: `P_RSR`, `P_CR`, `P_BCR` and `p_rsr`, `p_cr`, `p_bcr`.

## 2 Requirements

All methods are implemented in Python 3.12. The experiments have been run on a Linux server with a 3.7 GHz CPU and 128 GB of memory.

The required packages are:

- NetworkX 3.4.2
- SciPy 1.15.3

The environment can be created with:

```bash
conda create -n cnri python=3.12
conda activate cnri
pip install networkx==3.4.2 scipy==1.15.3
```

## 3 Datasets

The runner uses two comma-separated input files with the following formats.

### 3.1 Uncertain graph file

Each row represents an uncertain edge and contains the two endpoint IDs and its existence probability:

```text
source_vertex,target_vertex,existence_probability
```

| SOURCE VERTEX | TARGET VERTEX | EXISTENCE PROBABILITY |
| :---: | :---: | :---: |
| 0 | 1 | 0.17962666585239317 |
| 1 | 2 | 0.8835222053851841 |

### 3.2 Expected common-neighbor file

Each row contains an edge and the expected number of common neighbors between its endpoints:

```text
source_vertex,target_vertex,expected_number_of_common_neighbors
```

| SOURCE VERTEX | TARGET VERTEX | EXPECTED NUMBER OF COMMON NEIGHBORS |
| :---: | :---: | :---: |
| 0 | 1 | 0.0 |
| 1 | 3 | 0.23425672562750455 |

## 4 Usage

### 4.1 Run all methods

```bash
conda activate cnri
python run_all_algorithms.py
```

The methods are run in the following order:

```text
PIN, P+RSR, P+CR, P+BCR, ADR, TRPW, ABM
```

### 4.2 Run selected methods

Use `--algorithms` followed by one or more method names:

```bash
python run_all_algorithms.py --algorithms p_rsr p_cr p_bcr
```

Valid values are:

```text
pin p_rsr p_cr p_bcr adr trpw abm
```

### 4.3 Parameters

| PARAMETER | DEFAULT | DESCRIPTION |
| --- | --- | --- |
| `--graph` | `dataset/email-Enron/email-Enron-p.txt` | Path of the uncertain graph file. |
| `--cnexp` | `dataset/email-Enron/email-Enron_cnexp.dat` | Path of the expected common-neighbor file. |
| `--output-dir` | `output/` | Directory used to store the output possible worlds. |
| `--algorithms` | All methods | Methods to execute. |
| `--seed` | `42` | Random seed. |
| `--separator` | `,` | Input-file separator. |
| `--p-rsr-factor` | `1.0` | Sets the P+RSR iteration count to `round(|E| × factor)`. |
| `--p-cr-factor` | `1.0` | Sets the P+CR iteration count to `round(|E| × factor)`. |
| `--adr-factor` | `10.0` | Sets the ADR iteration count to `round(sum(p(e)) × factor)`. |
| `--trpw-factor` | `20.0` | Sets the TRPW iteration count to `round(sum(p(e)) × factor)`. |
| `--th-q` | `0.02` | Quality threshold `tau_q` of P+BCR. |
| `--th-c` | `0.9` | Confidence threshold `tau_c` of P+BCR. |


## 5 Output

All output possible worlds are stored in the `output/` directory. No spreadsheet file is produced.

| METHOD | OUTPUT FILE |
| --- | --- |
| PIN | `output/email_pin_output.dat` |
| P+RSR | `output/email_p_rsr_output.dat` |
| P+CR | `output/email_p_cr_output.dat` |
| P+BCR | `output/email_p_bcr_output.dat` |
| ADR | `output/email_adr_output.dat` |
| TRPW | `output/email_trpw_output.dat` |
| ABM | `output/email_abm_output.dat` |

Each output file represents one deterministic possible world. Each row contains one undirected edge:

```text
source_vertex,target_vertex
```

The console reports the number of edges, mean absolute discrepancy (MAD), and running time of each method:

```text
[1/7] PIN ...
  done: edges=..., mad=..., runtime=...s
```

MAD is the total discrepancy of a possible world divided by the number of uncertain edges. A smaller MAD indicates higher refinement quality.

Running the script again with the same output directory overwrites output files with the same names.

## 6 Project Structure

```text
.
├── algorithm/
│   ├── __init__.py
│   ├── our_algorithms.py
│   └── baselines.py
├── dataset/
│   └── email-Enron/
│       ├── email-Enron-p.txt
│       └── email-Enron_cnexp.dat
├── output/
│   ├── email_pin_output.dat
│   ├── email_p_rsr_output.dat
│   ├── email_p_cr_output.dat
│   ├── email_p_bcr_output.dat
│   ├── email_adr_output.dat
│   ├── email_trpw_output.dat
│   └── email_abm_output.dat
├── README.md
└── run_all_algorithms.py
```
