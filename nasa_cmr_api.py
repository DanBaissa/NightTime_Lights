import requests
from datetime import timedelta

def search_nasa_cmr(collection_id, start_date, end_date, bounding_box):
    cmr_search_url = "https://cmr.earthdata.nasa.gov/search/granules.json"
    bbox_str = f"{bounding_box[0]},{bounding_box[1]},{bounding_box[2]},{bounding_box[3]}"
    h5_links = []
    current_date = start_date
    while current_date <= end_date:
        target_date = current_date.strftime('%Y-%m-%d')
        params = {
            "collection_concept_id": collection_id,
            "temporal": target_date,
            "bounding_box": bbox_str,
            "page_size": 50
        }
        response = requests.get(cmr_search_url, params=params)
        if response.status_code == 200:
            granules = response.json()['feed']['entry']
            for granule in granules:
                granule_date = granule['time_start'].split('T')[0]
                if granule_date == target_date:
                    links = granule.get('links', [])
                    for link in links:
                        if 'href' in link and link['href'].startswith('https') and link['href'].endswith('.h5'):
                            h5_links.append(link['href'])
        current_date += timedelta(days=1)
    return h5_links