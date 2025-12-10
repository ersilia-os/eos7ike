# Entry-way rules (eNTRy) for gram-negative bacteria

The eNTRy rules are a set of three simple criteria for predicting compound accumulation in gram-negative bacteria.
These rules were derived from a large dataset of compound accumulation in Escherichia coli and have been shown
to be effective in guiding the design of compounds with improved permeability and accumulation in gram-negative bacteria.
The eNTRy rules state that compounds are more likely to accumulate in gram-negative bacteria if they possess the following
properties: (1) a low globularity (≤ 0.25), (2) a low number of rotatable bonds (≤ 5), and (3) the presence of a primary amine group.

This model was incorporated on 2025-12-04.


## Information
### Identifiers
- **Ersilia Identifier:** `eos7ike`
- **Slug:** `entry-rules`

### Domain
- **Task:** `Annotation`
- **Subtask:** `Activity prediction`
- **Biomedical Area:** `Antimicrobial resistance`
- **Target Organism:** `Escherichia coli`
- **Tags:** `Antimicrobial activity`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `3`
- **Output Consistency:** `Fixed`
- **Interpretation:** Three binary outputs indicating whether the compound meets each of the eNTRy rules, namely low globularity,
low number of rotatable bonds, and presence of a primary amine group. A value of 1 indicates that the compound
meets the respective criterion, while a value of 0 indicates that it does not.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| rb | integer | high | Low flexibility (rotatable bonds lower or equal than 5) |
| glob | integer | high | Low globularity (lower or equal than 0.25) |
| primary_amine | integer | high | Determines if a molecule has a primary amine |


### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos7ike.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos7ike.zip)

### Resource Consumption
- **Model Size (Mb):** `1`
- **Environment Size (Mb):** `623`


### References
- **Source Code**: [https://github.com/HergenrotherLab/entry-cli](https://github.com/HergenrotherLab/entry-cli)
- **Publication**: [https://www.nature.com/articles/nature22308](https://www.nature.com/articles/nature22308)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2017`
- **Ersilia Contributor:** [miquelduranfrigola](https://github.com/miquelduranfrigola)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [GPL-3.0-or-later](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos7ike
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos7ike
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
