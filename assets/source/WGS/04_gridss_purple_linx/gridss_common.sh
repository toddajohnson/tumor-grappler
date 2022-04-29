# Variables used in this file
#job_nodes=""
#def_slot=""
#s_vmem=""
#thread_num=""
#jvm_mem=""
#run_job_indexes=""
#base_log_filename=""
#logdir=""
#subject_id=""
#joint_sample_labels=""
#joint_bams_list=""
#gridss_version=""
#use_gridss_suffix=""
#steps=""
#base_shdir=""
#base_logdir=""
#JOINT_CALLING_INFO_FILE=""
#GRIDSS_JOINT_sh=""

gridss_script(){
	qsub -tc 5 \
		-t ${run_job_indexes} \
		-pe def_slot ${def_slot} \
		-l s_vmem="${s_vmem}G" \
		-N ${base_log_filename} \
		-o ${logdir} \
		-j y \
		-cwd \
		-v PGL_THREAD_NUM=${thread_num},PGL_JVM_MEM=${jvm_mem},STEPS=${steps},RUN_DIR=${run_dir},JOB_NODES=${job_nodes},GRIDSS_VERSION=${gridss_version},USE_GRIDSS_SUFFIX=${use_gridss_suffix},SUBJECT_ID=${subject_id},JOINT_SAMPLE_LABELS="${joint_sample_labels}",JOINT_BAMS="${joint_bams_list}" \
		${base_shdir}/${GRIDSS_JOINT_sh}
}

submit_gridss_script(){
	echo "Attempting to parse current job from ${JOINT_CALLING_INFO_FILE}"
	
	while IFS=$'\t' read -r subject_id subject_index sample_id_N joint_labels joint_bams; do
		if [[ "${subject_id}" != "subject.id" ]] ; then
			logdir="${base_logdir}/subject.`printf %02d ${subject_index}`"
			mkdir -p ${logdir}
	
			if [[ "${steps}" == "preprocess" ]]; then
				IFS=',' read -r -a labels_array <<< "$joint_labels"
				IFS=' ' read -r -a joint_bams_array <<< "$joint_bams"
		
				run_job_indexes="1"
				
				for index in "${!labels_array[@]}"
				do
					joint_sample_labels=${labels_array[index]}
					joint_bams_list=${joint_bams_array[index]}
					base_log_filename="preprocess.${joint_sample_labels}"
					
					echo "Running GRIDSS preprocessing on ${subject_id} for sample ${joint_sample_labels}."
	
					gridss_script
				done
			elif [[ "${steps}" == "all" ]] || [[ "${steps}" == "setupreference" ]] || [[ "${steps}" == "assemble" ]] || [[ "${steps}" == "call" ]] || [[ "${steps}" == "annotation" ]] || [[ "${steps}" == "recall" ]]; then
				joint_sample_labels=$(echo "${joint_labels}" | tr "," " ")
				joint_bams_list="${joint_bams}"
				
				if [[ "${steps}" == "all" ]]; then
					run_job_indexes="1"
					echo "Running all GRIDSS steps for ${subject_id}."
					base_log_filename="allsteps.${subject_id}"		
				elif [[ "${steps}" == "assemble" ]]; then
					if [[ $job_nodes > 1 ]]; then
						run_job_indexes="1-${job_nodes}"
						echo "Running GRIDSS assembly using ${job_nodes} jobs for ${subject_id}."
						base_log_filename="assembly.multijob.${subject_id}"		
					elif [[ $job_nodes == 1 ]]; then
						run_job_indexes="1"
						echo "Running GRIDSS single job assembly for ${subject_id}."
						base_log_filename="assembly.${subject_id}"	
					fi
				elif [[ "${steps}" == "call" ]]; then
					run_job_indexes="1"
					echo "Running GRIDSS SV calling for ${subject_id}."
					base_log_filename="call.${subject_id}"
				elif [[ "${steps}" == "recall" ]]; then
					run_job_indexes="1"
					echo "Running GRIDSS SV re-calling for ${subject_id}."
					base_log_filename="recall.${subject_id}"
				elif [[ "${steps}" == "annotation" ]]; then
					run_job_indexes="1"
					echo "Running GRIDSS kraken2 viral integration and RepeatMasker annotation for ${subject_id}."
					base_log_filename="annotate.${subject_id}"
				fi
	
				gridss_script
			fi
		fi
	done < ${JOINT_CALLING_INFO_FILE}
	
	echo "Finished running GRIDSS for ${subject_id}."
}


resubmit_gridss_script(){
	while IFS=$'\t' read -r subject_id subject_index sample_id_N joint_labels joint_bams job_index; do
		if [[ "${subject_id}" != "subject.id" ]] ; then
			logdir="${base_logdir}/subject.`printf %02d ${subject_index}`"
			mkdir -p ${logdir}
	
			if [[ "${steps}" == "assemble" ]] ; then
				joint_sample_labels=$(echo "${joint_labels}" | tr "," " ")
				joint_bams_list="${joint_bams}"
				
				run_job_indexes="${job_index}"
				echo "Re-running GRIDSS assembly using ${job_nodes} jobs for ${subject_id}."
				base_log_filename="assembly.multijob.${subject_id}"		

				gridss_script
			fi
		fi
	done < ${JOINT_CALLING_INFO_FILE}
	
	echo "Finished re-running GRIDSS for ${subject_id}."
}