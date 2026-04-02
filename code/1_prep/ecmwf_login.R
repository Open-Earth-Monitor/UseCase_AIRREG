options(keyring_backend="file")

user = "jheisig@uni-muenster.de"

# AMS User ID and API Key for CAMS access
uid_cams = "16538"
key_cams = "d2cfc6dc-6fc7-4785-a19d-40476c8f07f2"
#message("CAMS PW: theWeatherisn1ce!")


# CCS User ID and API Key for ERA5 access
uid_era = "149688"
#key_era = "3797cbc9-6304-47ac-a3ab-44197b08c873"
key_era = "0955e95e-6f24-41b5-9954-be2cea590031"
#message("ERA5 PW: theWeatherisn1ce!")


if(!("ecmwfr" %in% keyring::keyring_list()$keyring)){
  keyring::keyring_create("ecmwfr", password = "theWeatherisn1ce!")
}

if (keyring::keyring_is_locked("ecmwfr")){
  keyring::keyring_unlock("ecmwfr", password = "theWeatherisn1ce!")
}

# url: https://cds.climate.copernicus.eu/api
# key: 0955e95e-6f24-41b5-9954-be2cea590031
