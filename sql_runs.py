import os 
import time 
import duckdb
import main_sim as sim

curr_dir = os.getcwd().lower() 
files_folder_path = os.path.join(curr_dir, "Forward projection") 
filename_full_df = os.path.join(files_folder_path, "df_incid.parquet") 
duck_db_tmp_dir = os.path.join(curr_dir, "duck_db_temp")
output_file_path = os.path.join(files_folder_path, "df_incid_summed.parquet") 

sum_vacc_unvacc = False
output_peak_times = True

def output_peak_times_code(filename_full_df, output_folder_path):
    
    params = sim.init_params()
    scenarios = params["scenarios"]
    
    scenario_ids_sql = ", ".join(
        f"'{scenarios[key][4]}'"
        for key in scenarios.keys()
    )
    
    df_max = duckdb.sql(f"""
        WITH daily AS (
            SELECT
                simulation_index,
                scenario,
                run_nr,
                horizon,
                SUM(value) AS daily_value
            FROM read_parquet('{filename_full_df}')
            WHERE scenario IN ({scenario_ids_sql})
            GROUP BY
                simulation_index,
                scenario,
                run_nr,
                horizon
        )
        SELECT
            simulation_index,
            scenario,
            run_nr,
            arg_max(horizon, daily_value) AS peak_day
        FROM daily
        GROUP BY
            simulation_index,
            scenario,
            run_nr
    """).df()
    
    df_max.to_parquet(os.path.join(output_folder_path, "df_max_timing.parquet"), index = False)
    
def summ_vacc_unvacc_code():

    os.makedirs(duck_db_tmp_dir, exist_ok=True)
    
    con = duckdb.connect()
    
    con.execute("SET preserve_insertion_order=false")
    con.execute("SET threads=2")
    con.execute("SET memory_limit='10GB'")
    con.execute("SET max_temp_directory_size='50GB'")
    con.execute(f"SET temp_directory='{duck_db_tmp_dir}'")
    
    con.execute(
        f"""
        COPY (
            SELECT
                scenario,
                simulation_index,
                age_group,
                horizon,
                run_nr,
                SUM(value) AS value
            FROM read_parquet('{filename_full_df}')
            GROUP BY
                scenario,
                simulation_index,
                age_group,
                horizon,
                run_nr
        )
        TO '{output_file_path}' (FORMAT PARQUET)
        """
    )
    
    con.close()

if __name__ == "__main__":
    
    start_time = time.time() 
    
    if sum_vacc_unvacc:
        summ_vacc_unvacc_code()
    
    if output_peak_times:
        output_peak_times_code(filename_full_df, files_folder_path)
    
    end_time = time.time() 
    elapsed = end_time - start_time 
    print(f"Total time for running sql: {elapsed:.3f} seconds\n")


