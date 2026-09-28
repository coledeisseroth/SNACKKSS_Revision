#To run this, you must have a completed SNACKKSS_NLP pipeline run, located in your home directory. Feel free to replace the "~"s with wherever you keep that pipeline.
#Files needed for study and sample classification
for task in study sample; do
mkdir -p $task
for pert in gene drug; do
mkdir -p $task/$pert
for db in SNACKKSS_MC CREEDS; do
for stage in training testing; do
cp -r ~/SNACKKSS_NLP/$task/$pert/$db/biomedbert/${stage}_datasets/ $task/$pert/${db}_${stage}_datasets
done
done
done
done

#Target classification
mkdir -p target
for pert in gene drug; do
mkdir -p target/$pert
for db in CREEDS SNACKKSS_MC; do
cp -r ~/SNACKKSS_NLP/target/$pert/$db/biomedbert/training_datasets/ target/$pert/${db}_training_datasets
cp -r ~/SNACKKSS_NLP/target/$pert/$db/biomedbert/position_labels_test/ target/$pert/${db}_testing_datasets
done
done

#Control classification
mkdir -p control
for pert in gene drug; do
mkdir control/$pert
for db in CREEDS SNACKKSS_MC; do
cp -r ~/SNACKKSS_NLP/control/$pert/$db/biomedbert/training_json/ control/$pert/${db}_training_datasets
cp -r ~/SNACKKSS_NLP/control/$pert/$db/biomedbert/testing_datasets/ control/$pert/${db}_testing_datasets
done
done

#Gather the CREEDS-trained models
for stage in study sample target control; do
for pert in gene drug; do
cp -r ~/SNACKKSS_NLP/$stage/$pert/CREEDS/biomedbert/models/ $stage/$pert/CREEDS_models
done
done

