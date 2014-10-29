function sendEmail(attachment)
	if(isempty(ls('emailConfig.json')))
		disp('No emailConfig.json file found. Email will not be sent.');
		return;
	end

	emailConfig = loadjson('emailConfig.json');

	setpref('Internet','E_mail', emailConfig.senderEmailAddress);
	setpref('Internet','SMTP_Server','smtp.gmail.com');
	setpref('Internet','SMTP_Username', emailConfig.senderEmailAddress);
	setpref('Internet','SMTP_Password',emailConfig.senderEmailPassword);

	props = java.lang.System.getProperties;
	props.setProperty('mail.smtp.auth','true');
	props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
	props.setProperty('mail.smtp.socketFactory.port','465');

	tm = datestr(clock, 'YYYY-mm-dd hh:MM:ss');

	sendmail(emailConfig.recipientEmailAddress, ...
		'[SUCCESS] Turbulence_MonteCarlo Results', ['Email sent at ', tm], attachment);
end