sp_nm = all_sp_nm[1]
for (sp_nm in all_sp_nm){
  sp_nm = as.character(sp_nm) 
  sp_dir = paste0(str_replace_all(sp_nm,"_", "."), "/")
  temp_sp_files_to_delete<-paste0(project_path, sp_dir, "delete_temp_sp_files/", "*")
  unlink(temp_sp_files_to_delete, recursive = T, force = T) #delete previous frames
  # Loc <- "mydir"
  # system(paste0("rm -r ", temp_sp_files_to_delete))
}

temp_loc_to_delete = paste0(project_path, "temp/", "*")
unlink(temp_loc_to_delete, recursive = T, force = T) # delete previous frames

# system(paste0("rm -r ", temp_loc_to_delete))
