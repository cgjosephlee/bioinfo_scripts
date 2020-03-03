#!/usr/bin/env python3
'''
https://api.slack.com/apps
'''

import argparse
import sys
import os
import requests
import json

parser = argparse.ArgumentParser(description='Notify you in slack.')
parser.add_argument('msg', type=str,
                    help='message')
parser.add_argument('--user', type=str, default='hsinhan',
                    help='user name')
args = parser.parse_args()

user_name = args.user
msg = args.msg

token_fp = '/mnt/nas1/hhl/bin/My_scripts/slack_token.json'
with open(token_fp) as f:
    j = json.load(f)
    slack_token = j['token']

# get user ids
# param = {
#     'token': slack_token
# }
# r = requests.get('https://slack.com/api/users.list', params=param)
# with open('slack_users.json', 'w') as j:
#     json.dump(r.json(), j, indent=4)

user_json_fp = '/mnt/nas1/hhl/bin/My_scripts/slack_users.json'
user_id = None
with open(user_json_fp) as f:
    j = json.load(f)
    for item in j['members']:
        if item['profile']['display_name'] == user_name:
            user_id = item['id']
            break

if not user_id:
    raise OSError('User id for "{}" is not found.'.format(user_name))

print('Sending message to {} ({})...'.format(user_name, user_id), file=sys.stderr)

msg = '({}@{}) {}'.format(
    os.environ.get('USER'),
    os.uname().nodename,
    msg
)
param = {
    'token': slack_token,
    'channel': user_id,
    'text': msg
}
r = requests.get('https://slack.com/api/chat.postMessage', params=param)
j = r.json()
if not j['ok']:
    print(json.dumps(j, indent=2))
    sys.exit(1)
