import os, shutil
from flask import Flask, request, render_template, url_for, redirect, send_from_directory, flash, jsonify, send_file, Response
from forms import BamQuery_search, Retrieve_results
import secrets
import sys
import zipfile
import io
import pathlib
sys.path.append('../')
from BamQuery import running_for_web
from celery import Celery
import random, time

app = Flask(__name__)

app.config['SECRET_KEY'] = '226f3ae1820d1919837515c79dfe4a84'
app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'


celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)

@app.route("/search", methods=['POST', 'GET'])
def search():
	form = BamQuery_search()
	form_2 = Retrieve_results()
	image_file = url_for('static', filename='favicon.png')
	
	if form.validate_on_submit():
		name_exp = form.name_query.data
		parent_dir = os.path.join(app.root_path, 'static/temps')
		path_search = os.path.join(parent_dir, name_exp)
		path_input = os.path.join(path_search, 'Input')

		try:
			os.mkdir(path_search)
			os.mkdir(path_input)
			
			if form.peptides.data:
				path_saved, state = save_peptides(form.peptides.data, path_input)
				
				if state == 'Error':
					shutil.rmtree(path_search)
					task_ids[name_exp] = 'No task'
					flash(f'Sorry, peptides.tsv file contains more than 100 entries !', 'error')
				else:
					path_to_bam_directories = os.path.join(app.root_path, 'static/BAM_directories.tsv')
					shutil.copy2(path_to_bam_directories, path_input)
					task = running_BamQuery.apply_async(args = [path_input, name_exp, form.strandedness.data, form.th_out.data])
					flash(f'Please do not close this page, you must wait for the query to be completed !', 'success')
					return render_template('results.html', title ='Results', image_file = image_file, name_exp = name_exp, task_id = task.id)

		except FileExistsError:
			flash(f'Name of Query { name_exp } already in use. Please change the name of the query !', 'error')

	if form_2.validate_on_submit():
		name_exp = form_2.name_query_to_retrieve.data
		path_output = os.path.join(app.root_path, 'static/temps',name_exp,'output')
		filename = 'data.zip'
		try:
			shutil.make_archive(path_output, 'zip', path_output)
			path_output = path_output+'.zip'

			with open(path_output, 'rb') as f:
				data = f.readlines()
			shutil.rmtree(os.path.join(app.root_path, 'static/temps',name_exp))
			return Response(data, headers={
				'Content-Type': 'application/zip',
				'Content-Disposition': 'attachment; filename=%s;' % filename
			})
		except FileNotFoundError:
			flash(f'No query with this name: { name_exp } is found in our server !', 'error')	

		#return send_file(path_output, mimetype='application/zip', as_attachment=True, attachment_filename='data.zip')

	return render_template('search.html', title ='Search', image_file = image_file, form = form, form_2 = form_2, block=False)

def file_len(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

def save_peptides(form_peptides_file, path_input):
	peptides_file_path = os.path.join(path_input, form_peptides_file.filename)
	form_peptides_file.save(peptides_file_path)
	peptides_number = file_len(peptides_file_path)
	if peptides_number > 100:
		return peptides_file_path, 'Error'
	else:
		return peptides_file_path, 'Success'

@app.route('/')
@app.route('/about')
def about():
	return render_template('about.html')

@app.route('/documentation', methods=['POST', 'GET'])
def documentation():
	return redirect("http://bamquery.iric.ca/", code=302)

@app.route('/favicon.ico')
def favicon():
	return send_from_directory(os.path.join(app.root_path, 'static'), 'favicon.ico',mimetype='image/vnd.microsoft.icon')

@celery.task(bind=True)
def running_BamQuery(self, path_input, name_exp, strandedness, th_out):
	path_output = running_for_web(path_input, name_exp, strandedness, th_out)
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
	path_heatmap = os.path.join(app.root_path, 'static/temps',name_exp,'output','plots/heat_maps/r_plot',name_plot)
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
	shutil.rmtree(os.path.join(app.root_path, 'static/temps',name_exp))
	return redirect(url_for('search'))


if __name__ == '__main__':
	app.run(debug=True)
	port = int(os.environ.get('PORT', 33000))
	app.run(host='0.0.0.0', port=port)

