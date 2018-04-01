//
//  VitalSignsHandler.h
//  Pods
//
//  Created by Мария Тимофеева on 23.12.2017.
//

#import <Foundation/Foundation.h>

@interface VitalSignsHandler : NSObject
-(NSDictionary<NSString*, NSNumber*>*)addECG:(int)data;
-(NSDictionary<NSString*, NSNumber*>*)addPPGRed:(double)Red IR:(double)IR;
@end
