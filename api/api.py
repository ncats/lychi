from flask import Flask, request, jsonify, render_template
import subprocess
import tempfile
import os

lychify = Flask("lychify")

@lychify.route('/')
def index():
    return render_template('index.html')

@lychify.route('/standardize', methods=['POST'])
def standardize():
    # Get the list of SMILES objects from the request
    smiles_list = request.json['smiles_list']
    counter = 0
    if len(smiles_list) == 0:
        return jsonify({'result': []})
    # Write the SMILES objects to a temporary file
    with tempfile.NamedTemporaryFile(delete=False, mode='w') as f:
        for obj in smiles_list:
            if ("id" in obj):
                f.write(obj['smiles'] + '\t' + str(obj['id']) + '\n')
            else:
                counter += 1
                f.write(obj["smiles"] + '\t' + "no_id:" + str(counter) + '\n')
        filepath = f.name

    # Create the command to call the Java application with the input file
    command = ['java', '-jar', 'lychi-0.7.1-jar-with-dependencies.jar', filepath]

    # Call the Java application using subprocess
    process = subprocess.run(command, capture_output=True, text=True)

    # Remove the temporary file
    os.remove(filepath)

    # Check if the command completed successfully
    if process.returncode != 0:
        return jsonify({'error': 'Failed to run Java application',
                        'output': process.stdout.strip().split('\n')})

    # Parse the output from the Java application and format the response
    output = process.stdout.strip().split('\n')
    formatted_output = []
    for i, line in enumerate(output):
        smiles, id, lychi = line.split('\t')
        formatted_output.append({'smiles': smiles, 'id': id, 'lychi': lychi})

    # Return the formatted output as a JSON response
    return jsonify({'result': formatted_output})


if __name__ == '__main__':
    lychify.run(debug=True)
