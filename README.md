# tcrdist3_prep
Little script to prepare some of my TRUST4 outputs for use in tcrdist3. tcrdist3 seems to require the TCR tables squashed into one dataframe, so this sucks the files out by metadata, binds the patient ID as an extra column, concatenates all the tables and then renames the relevant fields to the names recognised by tcrdist3. First go using Spyder, which is a nice little environment for coding, has worked better for me than VS Code.
