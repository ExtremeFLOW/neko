
mypath=(""
"")
account=""

for i in ${!mypath[@]}; do
  input_folder_path=${mypath[$i]}
  case_num="${input_folder_path:74:3}"
  input_file_name="nusselt_t500.txt"
  input_file_path="${input_folder_path}/${input_file_name}"
  output_name="${case_num}-${input_file_name}"
  echo "secure copying $input_file_path"
  scp -r ${account}:${input_file_path} ${output_name}
done
