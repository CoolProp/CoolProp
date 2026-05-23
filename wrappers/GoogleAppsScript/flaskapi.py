from flask import Flask
from flask_restful import Api, Resource
import CoolProp as cp

app = Flask(__name__)
api = Api(app)

class PropsSI(Resource):
    def get(self, output, property_1, value_1, property_2, value_2,  fluid):
        output_value = cp.CoolProp.PropsSI(output, property_1, float(value_1), property_2, float(value_2), fluid)
        return output_value

class HAPropsSI(Resource):
    def get(self, output, property_1, value_1, property_2, value_2, property_3, value_3):
        output_value = cp.HumidAirProp.HAPropsSI(output, property_1, float(value_1), property_2, float(value_2), property_3, float(value_3))
        return output_value
    
@app.route('/')
def index():
    return 'CoolProp API'

api.add_resource(PropsSI, '/PropsSI/<output>/<property_1>/<value_1>/<property_2>/<value_2>/<fluid>')
api.add_resource(HAPropsSI, '/HAPropsSI/<output>/<property_1>/<value_1>/<property_2>/<value_2>/<property_3>/<value_3>')

if __name__ == '__main__':
    app.run()
    