# rxncon

The purpose of rxncon is to provide a framework to collect, visualise and model experimental data on cellular networks. In the rxncon framework, cellular signal transduction networks are described at the same granularity as empirical data. The key feature is strict separation of elemental reactions from contingencies, which define contextual constrains on these reactions, and this separation minimises the combinatorial complexity. The user defines the network as one reaction list and one contingency list. From these data mathematical and graphical representation can be generated. The network can be easily modified and extended, and both visualization and mathematical models can be generated automatically at any time.

For more details we refer to the following publications:

* Tiger, C.-F., Krause, F., Cedersund, G., Palmér, R., Klipp. E., Hohmann, S., Kitano, H. & Krantz, M. (2012)
A framework for mapping, visualisation and automatic model creation of signal transduction networks. [Molecular Systems Biology 8, 578](http://www.nature.com/msb/journal/v8/n1/full/msb201212.html).

* Romers, J.C. & Krantz, M. (2017) rxncon 2.0: a language for executable molecular systems biology. [bioRxiv:107136](https://doi.org/10.1101/107136)

## Installation

This software requires Python 3.5 or higher. Installation is straightforward
through the Python Package Index. 
If the `pip` command is linked to your Python 3 installation, run `pip install rxncon`,
otherwise run `pip3 install rxncon`.

This should install all libraries and the command-line tools to interface with them
on your machine.

## Usage

Please clone our [models repository](https://github.com/rxncon/models) to find
an example model describing the insulin pathway, as well as an example Excel
sheet that you can use to create your own models.

In what follows we assume the Excel file containing your network description
is called `model.xls`

### Compiling to a boolean network

To compile the rxncon system to a boolean network, run the command
```
rxncon2boolnet.py model.xls
``` 
This will generate three files `model.boolnet`, `model_symbols.csv` and `model_initial_vals.csv`

### Compiling to a rule-based model


### Generating graphs

## For devs

We accept pull requests. Some things to keep in mind:

* Please write tests for your code. The directory structure in `rxncon/test/` mirrors
the structure in the `rxncon/` directory, and we try to cover each module. Feel free
to look around. The testing framework we use is [pytest](https://docs.pytest.org/en/latest/)

* Please provide some documentation for each module, class, method and function. Don't overdo it either:
one-liner functions don't require 5 lines of documentation.

* Please use Python 3.5's [type annotations](https://www.python.org/dev/peps/pep-0484/) in your code. The shell script
`typecheck.sh` runs the [mypy](http://mypy-lang.org) type checker with the strictest options. This should give no errors.