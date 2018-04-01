//
//  VitalSignsHandler.m
//  Pods
//
//  Created by Мария Тимофеева on 23.12.2017.
//

#import "VitalSignsHandler.h"
#import "BridgeClass.h"
@interface VitalSignsHandler()
@property (strong, nonatomic) BridgeClass *bridge;
@end

@implementation VitalSignsHandler

- (instancetype)init
{
    self = [super init];
    if (self) {
        self.bridge = [BridgeClass new];
    }
    return self;
}

-(NSDictionary<NSString*, NSNumber*>*)addECG:(int)data{
  return [self.bridge updateECG:data];
}
-(NSDictionary<NSString*, NSNumber*>*)addPPGRed:(double)Red IR:(double)IR{
    return [self.bridge updatePPGRed:Red IR:IR];
}
@end
