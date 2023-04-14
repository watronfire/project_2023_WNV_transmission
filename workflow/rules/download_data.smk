rule download_sd_cases:
    message: "Downloads reported West Nile virus cases in California and San Diego. Data is reported https://data.sandiegocounty.gov/dataset/Department-of-Environmental-Health-Locally-Acquire/p9ku-563v"
    output:
        cases = "data/sd_cases.csv"
    run:
        cases = pd.read_csv( "https://data.sandiegocounty.gov/resource/p9ku-563v.csv", parse_dates=["date"] )
        cases.to_csv( output.cases, index=False )

rule combine_mosquito_data:
    message: "Combines mosquito data collected in three parts from San Diego Vector Control."
    input:
        batch1 = "data/temp1.csv",
        batch2 = "data/temp2.csv",
        batch3 = "data/temp3.csv"
    output:
        mosquito_numbers = "data/combined_mosquito_data.csv"
    shell:
        "python workflow/scripts/combine_mosquito_data.py"