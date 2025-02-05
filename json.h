/*
    Pheniqs : PHilology ENcoder wIth Quality Statistics
    Copyright (C) 2018  Lior Galanti
    NYU Center for Genetics and System Biology

    Author: Lior Galanti <lior.galanti@nyu.edu>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PHENIQS_JSON_H
#define PHENIQS_JSON_H

#include "include.h"
#include "error.h"
#include "kstring.h"

class ValidationError : public Error {
    public:
        Document ontology;
        ValidationError(Document& ontology) :
            Error("JSON directive validation error", ErrorCode::JSON_VALIDATION_ERROR),
            ontology(std::move(ontology)) {
            compile();
        };
        ValidationError(const ValidationError& other) :
            Error("JSON directive validation error", ErrorCode::JSON_VALIDATION_ERROR),
            ontology(kObjectType) {
            ontology.CopyFrom(other.ontology, ontology.GetAllocator());
            compile();
        };
        ostream& describe(ostream& o) const override;

    private:
        void compile();
};
Document encode_validation_error(const SchemaValidator& validator, const Value& schema, const Value& container);

void print_json(const Value& node, const char* path, const int32_t& precision=PHENIQS_FLOAT_PRECISION);
void print_json(const Value& node, ostream& o=cout, const int32_t& precision=PHENIQS_FLOAT_PRECISION);
// Document* load_json(const string& path);

/*  Recursively merge two JSON documents.
    Overlay < base > on < ontology >, using < document > for memory allocation.
    ontology is updated in the process.
*/
void merge_json_value(const Value& base, Value& ontology, Document& document);
void project_json_value(const Value& base, const Value& ontology, Value& container, Document& document);
void clean_json_value(Value& ontology, Document& document);
void sort_json_value(Value& ontology, Document& document);
void clean_json_object(Value& ontology, Document& document);
void overlay_json_object(Document& ontology, const Value& overlay);
bool remove_disabled_from_json_value(Value& ontology, Document& document);

/* decoding JSON container into an object */
template < typename T > bool decode_value(T& value, const Value& container);
template < typename T > bool decode_value_by_key(const Value::Ch* key, T& value, const Value& container);
template < typename T > T decode_value(const Value& container);
template < typename T > T decode_value_by_key(const Value::Ch* key, const Value& container);
// template < typename T > T decode_value(const Value& container) {
//     T value;
//     if(!decode_value(value, container)) {
//         throw ConfigurationError("value decoding failed");
//     }
//     return value;
// };
// template < typename T > T decode_value_by_key(const Value::Ch* key, const Value& container) {
//     T value;
//     if(container.IsObject()) {
//         if(!decode_value_by_key(key, value, container)) {
//             throw ConfigurationError("unknown key " + string(key));
//         }
//     } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
//     return value;
// };

/* encoding object to JSON container */
inline bool encode_key_value(const string& key, const bool& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), Value(value).Move(), document.GetAllocator());
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return true;
};
inline bool encode_key_value(const string& key, const uint8_t& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(value < numeric_limits< uint8_t >::max()) {
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), Value(static_cast< unsigned>(value)).Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const vector< uint8_t >& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            Value array(kArrayType);
            for(auto& v : value) {
                array.PushBack(Value(static_cast< unsigned >(v)).Move(), document.GetAllocator());
            }
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const list< uint8_t >& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            Value array(kArrayType);
            for(auto& v : value) {
                array.PushBack(Value(static_cast< unsigned>(v)).Move(), document.GetAllocator());
            }
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const int32_t& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(value < numeric_limits< int32_t >::max()) {
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), Value(static_cast< int >(value)).Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const vector< int32_t >& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            Value array(kArrayType);
            for(auto& v : value) {
                array.PushBack(Value(v).Move(), document.GetAllocator());
            }
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const uint32_t& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(value < numeric_limits< uint32_t >::max()) {
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), Value(static_cast< unsigned>(value)).Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const vector< uint32_t >& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            Value array(kArrayType);
            for(auto& v : value) {
                array.PushBack(Value(static_cast< unsigned>(v)).Move(), document.GetAllocator());
            }
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const int64_t& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(value < numeric_limits< int64_t >::max()) {
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), Value(value).Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const vector< int64_t >& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            Value array(kArrayType);
            for(auto& v : value) {
                array.PushBack(Value(v).Move(), document.GetAllocator());
            }
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const list< int64_t >& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            Value array(kArrayType);
            for(auto& v : value) {
                array.PushBack(Value(v).Move(), document.GetAllocator());
            }
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const uint64_t& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(value < numeric_limits< uint64_t >::max()) {
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), Value(value).Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const vector< uint64_t >& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            Value array(kArrayType);
            for(auto& v : value) {
                array.PushBack(Value(v).Move(), document.GetAllocator());
            }
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const list< uint64_t >& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            Value array(kArrayType);
            for(auto& v : value) {
                array.PushBack(Value(v).Move(), document.GetAllocator());
            }
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const double& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(value != numeric_limits< double >::infinity()) {
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), Value(value).Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const vector< double >& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            Value array(kArrayType);
            for(auto& v : value) {
                array.PushBack(Value(v).Move(), document.GetAllocator());
            }
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const list< double >& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            Value array(kArrayType);
            for(auto& v : value) {
                array.PushBack(Value(v).Move(), document.GetAllocator());
            }
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const string& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            container.AddMember (
                Value(key.c_str(), key.size(), document.GetAllocator()).Move(),
                Value(value.c_str(),value.length(), document.GetAllocator()).Move(),
                document.GetAllocator()
            );
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const vector< string >& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            Value array(kArrayType);
            for(auto& v : value) {
                array.PushBack(Value(v.c_str(), v.length(), document.GetAllocator()).Move(), document.GetAllocator());
            }
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const list< string >& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(!value.empty()) {
            Value array(kArrayType);
            for(auto& v : value) {
                array.PushBack(Value(v.c_str(), v.length(), document.GetAllocator()).Move(), document.GetAllocator());
            }
            container.AddMember(Value(key.c_str(), key.size(), document.GetAllocator()).Move(), array.Move(), document.GetAllocator());
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};
inline bool encode_key_value(const string& key, const kstring_t& value, Value& container, Document& document) {
    if(container.IsObject()) {
        container.RemoveMember(key.c_str());
        if(ks_not_empty(value)) {
            container.AddMember (
                Value(key.c_str(), key.size(), document.GetAllocator()).Move(),
                Value(value.s, value.l, document.GetAllocator()).Move(),
                document.GetAllocator()
            );
            return true;
        }
    } else { throw ConfigurationError(string(key) + " container is not a dictionary"); }
    return false;
};

#endif /* PHENIQS_JSON_H */
