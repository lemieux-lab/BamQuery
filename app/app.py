import os, shutil
from flask import Flask, request, render_template, url_for, redirect, send_from_directory, flash, jsonify, send_file, Response
from forms import BamQuery_search, Retrieve_results
import sys
sys.path.append('../')
from BamQuery import running_for_web
from celery import Celery
import random


app = Flask(__name__)

app.config['SECRET_KEY'] = '226f3ae1820d1919837515c79dfe4a84'
app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'


celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)

@app.route("/search", methods=['POST', 'GET'])
def search():
	form = BamQuery_search()
	image_file = url_for('static', filename='favicon.png')
	
	if form.validate_on_submit():
		name_exp = form.name_query.data.strip().replace(" ", "_")
		name_exp_aux = name_exp
		adding_random_name = str(random.getrandbits(16))
		name_exp = name_exp+'_'+adding_random_name

		parent_dir = os.path.join(app.root_path, 'static/temps')
		path_search = os.path.join(parent_dir, name_exp)
		path_input = os.path.join(path_search, 'Input')

		exists_search = os.path.exists(path_search) 
		exists_light = os.path.exists(path_input) 

		genome_version = form.genome_version.data
		data_base_snp = form.data_base_snp.data

		
		if form.peptides_list.data:
			if not exists_search and not exists_light:
				os.mkdir(path_search)
				os.mkdir(path_input)
			path_saved, state = save_peptides(form.peptides_list.data, name_exp_aux, path_input)
			path_to_bam_directories = os.path.join(app.root_path, 'static/BAM_directories.tsv')
			shutil.copy2(path_to_bam_directories, path_input)
			task = running_BamQuery.apply_async(args = [path_input, name_exp, True, form.th_out.data, genome_version, data_base_snp])
			message = 'Your query '+name_exp+' has been launched. \n Please do not close this page until your query has completed, otherwise you should get the query name to download the results later.'
			flash(message, 'success')
			return render_template('results.html', title ='Results', image_file = image_file, name_exp = name_exp, task_id = task.id)


	return render_template('search.html', title ='Search', image_file = image_file, form = form, block=False)

@app.route("/retrival", methods=['POST', 'GET'])
def retrival():
	form_2 = Retrieve_results()
	image_file = url_for('static', filename='favicon.png')
	
	if form_2.validate_on_submit():
		try:
			name_exp = form_2.name_query_to_retrieve.data
			path_output = os.path.join(app.root_path, 'static/temps',name_exp,'output')
			exists_output = os.path.exists(path_output)
			if exists_output:
				exists = os.path.exists(path_output+'/README.txt')
				if exists:
					filename = 'data.zip'
				
					shutil.make_archive(path_output, 'zip', path_output)
					path_output = path_output+'.zip'

					with open(path_output, 'rb') as f:
						data = f.readlines()
					try:
						os.system('rm -rf "{}"'.format(os.path.join(app.root_path, 'static/temps',name_exp, '/logs')))
						shutil.rmtree(os.path.join(app.root_path, 'static/temps',name_exp))
					except:
						pass
					return Response(data, headers={
						'Content-Type': 'application/zip',
						'Content-Disposition': 'attachment; filename=%s;' % filename
					})
				else:
					flash(f'The query { name_exp } is still in process !', 'error')
			else:
				flash(f'No query with this name: { name_exp } is found in our server !', 'error')
		except FileNotFoundError:
			flash(f'No query with this name: { name_exp } is found in our server !', 'error')

	return render_template('retrival.html', title ='Retrival', image_file = image_file, form_2 = form_2, block=False)


def file_len(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

def save_peptides(form_peptides_list, name_exp, path_input):
	to_write = ''
	peptides_list = form_peptides_list.split('\n')
	for pep in peptides_list:
		to_write+=pep.strip()+'\t'+name_exp+'\n'
	
	peptides_file_path = os.path.join(path_input, 'peptides.tsv')
	peptides_file = open(peptides_file_path, 'w')
	peptides_file.write(to_write)
	peptides_file.close()
	return peptides_file_path, 'Success'

@app.route('/')
@app.route('/about')
def about():
	return render_template('about.html')

@app.route('/installation')
def installation():
	return render_template('installation.html')

@app.route('/documentation', methods=['POST', 'GET'])
def documentation():
	return redirect("http://bamquery.iric.ca/", code=302)

@app.route('/favicon.ico')
def favicon():
	return send_from_directory(os.path.join(app.root_path, 'static'), 'favicon.ico',mimetype='image/vnd.microsoft.icon')

@celery.task(bind=True)
def running_BamQuery(self, path_input, name_exp, strandedness, th_out, genome_version, db_SNP):
	path_output = running_for_web(path_input, name_exp, strandedness, genome_version, db_SNP, th_out)
	self.update_state(state='PROGRESS', meta={'Query': name_exp, 'status': 'Running!'})
	return {'status': 'Finished!', 'result': 'Job completed!'}

@app.route('/status', methods=['POST'])
def status():
	name_exp = request.form['name_exp']
	task_id = request.form['task_id']
	return jsonify({}), 202, {'Location': url_for('taskstatus', name_exp = name_exp, task_id = task_id)}

@app.route('/status/<name_exp>/<task_id>')
def taskstatus(name_exp, task_id):
	task = running_BamQuery.AsyncResult(task_id)

	verb = ['Starting up', 'Booting', 'Repairing', 'Loading', 'Checking']
	adjective = ['master', 'radiant', 'silent', 'harmonic', 'fast']
	noun = ['solar array', 'particle reshaper', 'cosmic ray', 'orbiter', 'bit']
	message = '{0} {1} {2}...'.format(random.choice(verb), random.choice(adjective), random.choice(noun))

	if task.state == 'PENDING':
		# job did not start yet
		response = {
			'message': message,
			'Query': name_exp,
			'status': 'Running...',
			'state': task.state
		}
	elif task.state != 'FAILURE':
		response = {
			'message': message,
			'Query': task.info.get('Query', ''),
			'status': task.info.get('status', ''),
			'state': task.state
		}
		if 'result' in task.info:
			response['result'] = task.info['result']
	else:
		# something went wrong in the background job
		response = {
			'message': message,
			'Query': name_exp,
			'status': str(task.info),  # this is the exception raised
			'state': task.state
		}
	return jsonify(response)

@app.route("/display_heatmap/<name_exp>", methods=["GET"], endpoint='display_heatmap')
def display_heatmap(name_exp) :
	name_plot = name_exp+'_rna_norm_all_tissues.pdf'
	path_heatmap = os.path.join(app.root_path, 'static/temps',name_exp,'output','plots/heat_maps/transcription_evidence_heatmap/average_transcription_expression_heatmap',name_plot)
	return send_file(path_heatmap, attachment_filename = name_plot)

@app.route("/display_biotype/<name_exp>", methods=["GET"], endpoint='display_biotype')
def display_biotype(name_exp) :
	name_plot = name_exp+'_All_peptides.pdf'
	path_heatmap = os.path.join(app.root_path, 'static/temps',name_exp,'output','plots/biotypes/biotype_by_sample_group/all_peptides',name_plot)
	return send_file(path_heatmap, attachment_filename = name_plot)

@app.route('/download-zip/<name_exp>', methods=["GET"], endpoint='download-zip')
def request_zip(name_exp):
	path_output = os.path.join(app.root_path, 'static/temps',name_exp,'output')
	shutil.make_archive(path_output, 'zip', path_output)
	path_output = path_output+'.zip'
	return send_file(path_output, mimetype='application/zip', as_attachment=True, attachment_filename='data.zip')

@app.route('/remove', methods=['POST'])
def remove():
	name_exp = request.form['name_exp']
	try:
		shutil.rmtree(os.path.join(app.root_path, 'static/temps',name_exp))
	except:
		pass
	return redirect(url_for('search'))


if __name__ == '__main__':
	port = int(os.environ.get('PORT', 35000))
	app.run(host='0.0.0.0', port=port, debug=True)

