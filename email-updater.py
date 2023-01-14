### Jonathan Curtis
### 01/14/2022
### Test email updater

import smtplib, ssl

smtp_server = "smtp.gmail.com"
port = 465

sender_email= ""###Your email here
password = input("Type password: ")

receiver_email = ""###Phone number @ portal here

message = """Subject: Calculation

Calculation n/N done, runtime: s"""

context = ssl.create_default_context()

with smtplib.SMTP_SSL(smtp_server,port,context=context) as server:
	server.login(sender_email,password)
	server.sendmail(sender_email,receiver_email,message)