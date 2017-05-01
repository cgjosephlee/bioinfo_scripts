#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Im_Uaena.py
Get image urls from tistory fancam blogs, and upload to imgur

author: JLee
date: 2017.4
"""

from __future__ import print_function
import sys
import os
import argparse
import re
import requests
from bs4 import BeautifulSoup
# from selenium import webdriver # dummy browser
from imgurpython import ImgurClient
from imgurpython.helpers.error import ImgurClientError
from imgurpython.helpers.error import ImgurClientRateLimitError
import json
# from pprint import pprint

# test arguments
#sys.argv = ["prog", "-n", "97", "170402.url.txt", "test album"]
#sys.argv = ["prog", "-n"]

# client info
user_name = "cgjosephlee"
token_path = "~/" + user_name + ".token.json"
blacklist_path = "~/blacklist.txt"


def parse_arg():
    parser = argparse.ArgumentParser(
            description='Get image urls from tistory fancam blogs, and '
                        'upload to imgur.')
#            usage='{} url_list [album_title/url]'.format(sys.argv[0]))
    parser.add_argument('url_list', type=unicode,
                        help='List of tistory blog urls')
    parser.add_argument('album', type=unicode, metavar='album_title/url',
                        nargs='?', default='',
                        help='Can be either an album title '
                             'or an existing ablum url, upload as non-album '
                             'image if not spedifed')
    parser.add_argument('-p', required=False, action='store_true',
                        help='Parse image urls only, not upload')
    parser.add_argument('-l', required=False, action='store_true',
                        help='Assuming urls in the list are direct '
                             'image links')
    parser.add_argument('-u', required=False, action='store_true',
                        help='Print imgur urls when finish uploading')
    parser.add_argument('-n', metavar="NUM", default='0', type=int,
                        help='Ignore several images from head')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    # return vars(parser.parse_args())  # save as dict
    return parser.parse_args()


def authentication():
    global client, acc_token
    # open token file
    try:
        with open(os.path.expanduser(token_path)) as f:
            data = json.load(f)
            client_id = data['client_id']
            client_secret = data['client_secret']
            acc_token = data['access_token']
            ref_token = data['refresh_token']
        client = ImgurClient(client_id, client_secret,
                             acc_token, ref_token)
    except IOError:
        print("Error: token file not find")
        sys.exit(1)
    except ImgurClientError as e:
        print(e.status_code, e.error_message)
        sys.exit(1)


def write_auth():
    if acc_token != client.auth.current_access_token:
        token_json = {
                'user_name': user_name,
                'client_id': client.client_id,
                'client_secret': client.client_secret,
                'access_token': client.auth.current_access_token,
                'refresh_token': client.auth.refresh_token
                }
        with open(os.path.expanduser(token_path), "w") as f:
            json.dump(token_json, f, indent=4)


def get_img_link(url_list):
    headers = {"User-Agent": "curl"}
    re_tistory = r"http://cfile\d{1,2}.uf.tistory.com/image/(\w*)"
    img_url = []
    for url in url_list:
        web = requests.get(url, headers=headers)
        soup = BeautifulSoup(web.text, "html5lib")
        for link in soup.find_all("img"):
            raw_img = link.get("src")
            if re.match(re_tistory, raw_img) and \
               re.match(re_tistory, raw_img).group(1) not in blacklist:
                img_url.append(re.sub("image", "original", raw_img))
    return img_url


def get_album_id(arg):
    if arg == '':
        id = ''
    # uplaod to an exist album with specified album url
    # e.g. 'http://imgur.com/a/sOaVB'
    elif re.match(r'http://imgur.com/a/(\w{5})$', arg):
        id = re.match(r'http://imgur.com/a/(\w{5})$', arg).group(1)
    # given an album title
    else:
        album_list = client.get_account_albums('me')  # cost 1 credict
        album_titles = []
        for n in range(10):  # check 10 newest album titles only
            album_titles.append(album_list[n].title)
        if arg in album_titles:
            id = album_list[album_titles.index(arg)].id
        else:
            fields = {'title': arg, 'privacy': 'hidden'}
            new_album = client.create_album(fields)
            id = new_album['id']
    return id


def upload_img(img_list, album_id, n):
    if re.match(r'^\w{5}$', album_id):
        config = {'album': album_id}
    else:
        config = {}
    uploaded_img = []
    limit = False
    for url in img_list[n:]:
        for attempt in range(3):
            try:
                img = client.upload_from_url(url, config=config, anon=False)
                # cost 10 credicts
                uploaded_img.append(img['link'])
                # print('\r{}/{}'.format(n, len(img_list)), end='\r')
                sys.stdout.write('\r{}/{}'.format(n+1, len(img_list)))
                sys.stdout.flush()
                n += 1
            except ImgurClientError as e:
                print('\nError {}: {}. (Attempt {})'.format(e.status_code,
                      e.error_message, attempt+1))
            except ImgurClientRateLimitError:
                # to show the remaining waitng time of current IP
                r = requests.post("https://api.imgur.com/3/image.json",
                                  headers={'Authorization': 'Bearer ' +
                                           client.auth.current_access_token},
                                  data={'album': album_id,
                                        'image': url})
                print("\n"+r.json()['data']['error']['message'])
                limit = True
                break
            else:  # try clause completes normally
                break
        else:  # executes when the loop attempts 3 times
            print("\nAttempts is exceeded, maybe try later. ({}/{})".format(n,
                  len(img_list)))
            limit = True
        if limit:
            break
    else:  # executes when the loop completes normally
        print("\nUpload finished!")
    return uploaded_img


# def uniqList(seq):
#     # http://stackoverflow.com/a/480227
#     seen = set()
#     seen_add = seen.add
#     return [x for x in seq if not (x in seen or seen_add(x))]

# =============================================================================
args = parse_arg()

# open url list
try:
    with open(args.url_list) as f:
        url_list = []
        for line in f:
            line = line.strip()
            url_list.append(line)
except IOError:
    print("Error: can\'t find or read url list")
    sys.exit(1)

# read pre-defined blacklist
blacklist = []
try:
    with open(os.path.expanduser(blacklist_path)) as f:
        for line in f:
            line = line.strip()
            blacklist.append(line)
except IOError:
    print("No blacklist was imported")

# main
if args.l:
    img_url = url_list
else:
    img_url = get_img_link(url_list)

if args.p:
    for i in img_url:
        print(i)
else:
    authentication()
    album_id = get_album_id(args.album)
    uploaded_img = upload_img(img_url, album_id, args.n)
    if args.u and len(upload_img) > 0:
        for i in uploaded_img:
            print(i)
    write_auth()
