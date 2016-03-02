#!/bin/bash

# set -e

user_name="cgjosephlee"
client_id="516a14230cc31d8"
client_secret="5bde3ac197c579bb8fc3fadc54e7466d8cad05c3"

function usage {
    echo "Usage: $(basename $0) [-hn] list [album_title]" >&2
    echo "  -h  show this message"
    echo "  -n  urls in the list are direct image links"
    echo "  album_title can be either a new album title or an existing ablum url" >&2
    echo "  If no given album_title, it will output url directly" >&2
    exit 0
}

# defalut value
parsing=true

while getopts ":hn" options; do
    case $options in
        h)
        usage
        ;;
        n)
        parsing=false
        ;;
        \?)
        usage
        ;;
    esac
done

shift $((OPTIND -1))
list=$1
album_title=$2

if [ $# == 0 ]; then
    usage
elif [ ! -s $list ]; then
    echo "No list sepcified!"
    echo ""
    usage
fi

if $parsing; then
    echo "Fetching url..." >&2
    # download web source code from the list
    web_content=$(wget -q -O - -i $list)

    # fetch the image url from the source code
    img_url=$(echo $web_content | grep -oh 'http://cfile\d*.uf.tistory.com/image/\w*' | uniq | sed 's/image/original/g')
else
    echo "Assume urls in list are image links, no parsing performed..." >&2
    img_url=$(cat $list)
fi

num_img=$(echo $img_url | tr " " "\n" | wc -l)
echo "Found $num_img images..." >&2

# check album_title
re="http://imgur.com/a/([a-zA-Z0-9]{5})"
if [ ! "$album_title" ]; then
    echo "No album name sepcified, output url directly..." >&2
    echo $img_url | tr " " "\n"
    exit 0
elif [[ $album_title =~ $re ]]; then
    create_album=false
    album_id=${BASH_REMATCH[1]}
else
    create_album=true
fi

# get old token
if [ ! -s ${user_name}.token ]; then
    echo "Token file not found or broken!" >&2
    exit 1
else
    token=$(tail -n1 ${user_name}.token | sed -E 's/.*<access_token>(.*)<\/access_token>.*/\1/')
fi

# check if token expired
if [ -n $token ]; then
    # use "get settings" function to check if token expired
    response=$(curl -X GET -H "Authorization: Bearer $token" https://api.imgur.com/3/account/me/settings.xml 2> /dev/null)

    if [ $(echo $response | sed -E 's/.*status="(.*)".*/\1/g') == 403 ]; then # old token expired
        token_expire=true
        echo "Token expired, need to refresh token" >&2
    elif [ $(echo $response | sed -E 's/.*status="(.*)".*/\1/g') == 200 ]; then # success
        token_expire=false
    else # other error code
        err_code=$(echo $response | sed -E 's/.*status="(.*)".*/\1/g')
        err_msg=$(echo $response | sed -E 's/.*<error>(.*)<\/error>.*/\1/')
        echo "Error:" $err_code $err_msg >&2
        exit 1
    fi
fi

# regenerate token if old one expired
if $token_expire; then
    echo "Refreshing token..." >&2
    # get old refresh token
    refresh_token=$(tail -n1 ${user_name}.token | sed -E 's/.*<refresh_token>(.*)<\/refresh_token>.*/\1/')
    # gain new token and write to file
    curl -X POST -F "refresh_token=$refresh_token" -F "client_id=$client_id" -F "client_secret=$client_secret" -F "grant_type=refresh_token" https://api.imgur.com/oauth2/token.xml > ${user_name}.token 2> /dev/null
    # get new token from token file
    token=$(tail -n1 ${user_name}.token | sed -E 's/.*<access_token>(.*)<\/access_token>.*/\1/')
fi

# create an album
if $create_album; then
    echo "Creating new album..." >&2
    response=$(curl -X POST -H "Authorization: Bearer $token" -F "title=$album_title" https://api.imgur.com/3/album.xml 2> /dev/null)
    album_id=$(echo $response | sed -E 's/.*<id>(.*)<\/id>.*/\1/')
else
    echo "Upload to existing album..." >&2
fi

# upload image urls to the album
echo "Uploading images from urls..." >&2
for i in $img_url; do
    response=$(curl -X POST -F "album=$album_id" -F "image=$i" -H "Authorization: Bearer $token" https://api.imgur.com/3/image.xml 2> /dev/null)
    
    if [ $(echo $response | sed -E 's/.*status="(.*)".*/\1/g') == 200 ]; then # success
        link=$(echo $response | sed -E 's/.*<link>(.*)<\/link>.*/\1/')
        echo $link
    else # other error code
        err_code=$(echo $response | sed -E 's/.*status="(.*)".*/\1/g')
        err_msg=$(echo $response | sed -E 's/.*<error>(.*)<\/error>.*/\1/')
        echo "Error:" $err_code $err_msg >&2
        echo $i >&2
    fi
done
