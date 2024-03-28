from urllib.parse import urlencode


def get_download_url(ra:float, dec:float, pixel_size:int, bands:str='grz', data_release:int=10):
    """
    Generate the DECALS cutout service URL used to download each image, as either rgb or jpeg.

    Download urls with `cat urls.txt | shuf | xargs -n10 -P4 wget --continue` or similar. I find this quicker than all within Python.

    Args:
        ra (float): right ascension of galaxy
        dec (float): declination of galaxy
        pixel_size (int) pixels to download in native resolution
        bands (str): of bands to download e.g. grz, griz, gz, etc. Default is grz.
        data_release (str): DECALS data release to source image from
    Returns:
        (str): url to download galaxy in requested size/format
    """
    params = {
        'ra': ra,
        'dec': dec,
        'size': pixel_size,
        'bands': bands,
        'layer': f'ls-dr{data_release}'
    }
    assert isinstance(data_release, int)
    url = "http://legacysurvey.org/viewer/fits-cutout?" + urlencode(params)

    return url


if __name__ == '__main__':
    url = get_download_url(
        ra=150.0,
        dec=2.0,
        pixel_size=424,
        data_release=10,
    )
    print(url)
    # e.g. http://legacysurvey.org/viewer/fits-cutout?ra=150.0&dec=2.0&size=424&bands=grz&layer=ls-dr10
