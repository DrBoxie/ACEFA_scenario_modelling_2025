import os 
import time 
import duckdb 

start_time = time.time() 

curr_dir = os.getcwd().lower() 
files_folder_path = os.path.join(curr_dir, "Forward projection") 
filename_full_df = os.path.join(files_folder_path, "df_incid.parquet") 
output_file_path = os.path.join(files_folder_path, "df_incid_summed.parquet") 

duckdb.sql(
    f""" 
    COPY ( 
        SELECT 
        scenario, 
        simulation_index, 
        age_group, 
        horizon, 
        run_nr, 
        SUM(value) AS value 
        FROM '{filename_full_df}' 
        GROUP BY 
        scenario, 
        simulation_index, 
        age_group, 
        horizon, run_nr ) 
    TO '{output_file_path}' (FORMAT PARQUET) 
    """) 
    
end_time = time.time() 
elapsed = end_time - start_time 
print(f"Total time for summing vacc and unvacc: {elapsed:.3f} seconds\n")
