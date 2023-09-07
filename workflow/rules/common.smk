from snakemake.utils import validate
import pandas as pd

# with --use-conda --use-singularity

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

cancers = samples['cancer'].values.tolist()

cancersUniq = set(cancers)
cancersUniq = list(cancersUniq)

out_dir = config["out_dir"]
image_base_url = config["image_base_url"]
genome = config["genome"]
method = config["method"]
