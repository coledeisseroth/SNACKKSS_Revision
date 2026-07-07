for dir in cleaned study_comparison sample_comparison target_comparison control_comparison; do
cd $dir
bash BUILD.sh
cd ..
done

