import nimare
import pandas as pd
from nimare.dataset import Dataset

# First, load the dataset as `dset`

# Read in the topic file, rename the ID column, and
# prepend a prefix to the topic names
df = pd.read_table("/data/hu_fan/Desktop/test/v3/v3-topics-50.txt")
topic_names = [c for c in df.columns if c.startswith("topic")]
topics_renamed = {t: "Neurosynth_LDA__" + t for t in topic_names}
topics_renamed["id"] = "study_id"
df = df.rename(columns=topics_renamed)

# Change the data type for the study_id column so it can be merged
df['study_id'] = df['study_id'].astype(str)

dset = Dataset.load("/data/hu_fan/Desktop/test/neurosynth_dataset_with_abstracts.pkl.gz")


# Merge the topic dataframe into the annotations dataframe
new_annotations = dset.annotations.merge(
    df, 
    how="inner", 
    left_on="study_id", 
    right_on="study_id"
)
#change dset.annotations
dset.annotations = new_annotations

# The topic file only contains ~10k studies,
# so we must reduce the dataset to match
new_ids = new_annotations["id"].tolist()
# do overall reductions--9950
dset = dset.slice(new_ids)
