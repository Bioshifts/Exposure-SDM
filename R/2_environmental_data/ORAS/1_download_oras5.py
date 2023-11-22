import cdsapi

c = cdsapi.Client()

for year in range(1958, 2020):
    c.retrieve(
        'reanalysis-oras5',
        {
            'product_type': [
                'consolidated', 'operational',
            ],
            'vertical_resolution': 'single_level',
            'variable': 'sea_surface_temperature',
            'year': 
                str(year),
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'format': 'zip',
        },
        '/scratch/boliveira/Data/Marine/oras/SST_ORAS_downloads/download_'+str(year)+'.zip')
                
