import os
from pprint import pprint

import nimare

import nibabel as nib
import numpy as np
from nilearn import image, masking, plotting

from nimare import annotate, decode
from nimare.dataset import Dataset
from nimare.utils import get_resource_path


from nimare.extract import download_abstracts, fetch_neuroquery, fetch_neurosynth
from nimare.io import convert_neurosynth_to_dataset

#### generating topic-based Dataset

out_dir = os.path.abspath("/data/p_02542/Fan_ThalamicProject/neurosynth/CPcTFCE/test/")
os.makedirs(out_dir, exist_ok=True)

# Fetch Neurosynth with *just* the LDA50 features
files = nimare.extract.fetch_neurosynth(
    data_dir=out_dir,  # version 0.0.10 switched to data directory
    version="7",
    overwrite=False,
    source="abstract",
    vocab="LDA50",  # Note the difference here
)
neurosynth_db = files[0]
pprint(neurosynth_db)
# Note the "keys" file. That has the top 30 words for each topic.
# It *doesn't* go in the Dataset at all though.

# Get the Dataset object
neurosynth_dset = nimare.io.convert_neurosynth_to_dataset(
    coordinates_file=neurosynth_db["coordinates"],
    metadata_file=neurosynth_db["metadata"],
    annotations_files=neurosynth_db["features"],
)
neurosynth_dset.save('Dataset.pkl')

# Add article abstracts to dataset
neurosynth_dset = download_abstracts(neurosynth_dset, "fan@cbs.mpg.de")
neurosynth_dset.save(os.path.join(out_dir, "neurosynth_dataset_with_abstracts.pkl.gz"))







