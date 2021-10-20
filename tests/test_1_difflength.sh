echo diff_length linear model, expect table of results with mean_diff
../diff_length.R -d test_data/plen_test_data.tab -b "control" -m test_data/plen_test_metadata.tab

echo diff_length log linear model, expect table of results with log2FC
../diff_length.R -d test_data/plen_test_data.tab -b "control" -l TRUE -m test_data/plen_test_metadata.tab

echo diff_length mixed model, expect table of results with mean_diff
../diff_length.R -d test_data/plen_test_data.tab -p zee -b "control" -l FALSE -t m -m test_data/plen_test_metadata.tab

echo diff_length mixed model, expect table of results with log2FC
../diff_length.R -d test_data/plen_test_data.tab -p zee -b "control" -l TRUE -t m -m test_data/plen_test_metadata.tab

echo diff_length using missing variables, expect stop
../diff_length.R -d test_data/plen_test_data.tab -p vee+tee*dee -m test_data/plen_test_metadata.tab

echo diff_length using Wilcoxon, expect table of results for Wilcoxon
../diff_length.R -d test_data/plen_test_data.tab -b "control" -t w -m test_data/plen_test_metadata.tab

echo diff_length using Wilcoxon, expect additional warning
../diff_length.R -d test_data/plen_test_data.tab -p zee -b "control" -t w -m test_data/plen_test_metadata.tab
